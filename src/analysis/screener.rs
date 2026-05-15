//! Screening logic using pairwise alignment.
//!
//! Iterates through the template sequence with different oligo lengths,
//! using pairwise alignment to find best matches in each reference sequence.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;

use rayon::prelude::*;

use super::analyzer::analyze_sequences;
#[cfg(feature = "cuda")]
use super::cuda_align::{
    CudaAligner, collect_matches_with_cuda_aligner,
    collect_mismatch_counts_with_cuda_aligner, create_cuda_aligner,
};
use super::fasta::{ReferenceData, TemplateData};
use super::pairwise::{
    DnaAligner, collect_matches_with_aligner, collect_mismatch_counts_with_aligner, create_aligner,
};
use super::simd_align::{
    SimdAligner, collect_matches_with_simd_aligner,
    collect_mismatch_counts_with_simd_aligner, create_simd_aligner,
};
use super::simple_align::{
    SimpleAligner, collect_matches_with_simple_aligner,
    collect_mismatch_counts_with_simple_aligner, create_simple_aligner,
};
use super::types::{
    AlignerKind, AnalysisParams, ExclusivityResult, LengthResult, MismatchBucket, PositionResult,
    ProgressUpdate, ScreeningResults, WindowAnalysisResult,
};
use crate::pause::PauseFlag;

/// Per-worker aligner state, holding whichever backend is in use.
enum WorkerAligner {
    Pairwise(DnaAligner),
    Simple(SimpleAligner),
    SimpleSimd(SimdAligner),
    #[cfg(feature = "cuda")]
    SimpleCuda(CudaAligner),
}

/// Build the list of oligo lengths to process given min, max, and skip.
///
/// `skip = 0` produces every length. `skip = N` skips N lengths between
/// processed ones (so the step is `skip + 1`). `max` is included if it lands
/// on a step boundary.
fn build_length_list(min: u32, max: u32, skip: u32) -> Vec<u32> {
    if min > max {
        return Vec::new();
    }
    let step = skip.saturating_add(1).max(1) as usize;
    (min..=max).step_by(step).collect()
}

/// Run the complete screening analysis using pairwise alignment.
pub fn run_screening(
    template: &TemplateData,
    references: &ReferenceData,
    params: &AnalysisParams,
    exclusivity: Option<&ReferenceData>,
    progress_tx: Option<Sender<ProgressUpdate>>,
    pause: Option<PauseFlag>,
) -> ScreeningResults {
    let num_threads = params.thread_count.get_count();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

    let differential_enabled = exclusivity.is_some();
    let exclusivity_sequence_count = exclusivity.map(|e| e.len());

    let mut results = ScreeningResults::new(
        params.clone(),
        template.sequence.len(),
        references.len(),
        template.sequence.clone(),
        differential_enabled,
        exclusivity_sequence_count,
    );

    let ref_bytes: Vec<Vec<u8>> = references
        .sequences
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let ref_bytes = Arc::new(ref_bytes);

    let excl_bytes: Option<Arc<Vec<Vec<u8>>>> = exclusivity.map(|e| {
        Arc::new(e.sequences.iter().map(|s| s.as_bytes().to_vec()).collect())
    });
    let excl_names: Option<Arc<Vec<String>>> = exclusivity.map(|e| Arc::new(e.names.clone()));

    // If the GPU backend is selected, upload the reference sets once before
    // the per-position loop. Slot 0 is the main set; slot 1 is the exclusivity
    // set when differential mode is enabled. Per-position dispatches identify
    // their slot by slice pointer + length, so the same `Arc<Vec<Vec<u8>>>`
    // must be the one handed to each worker — the rest of this function
    // already does that.
    #[cfg(feature = "cuda")]
    if params.aligner == AlignerKind::SimpleCuda {
        use super::cuda_align::{ensure_initialized, register_slot};
        ensure_initialized()
            .unwrap_or_else(|e| panic!("CUDA backend init failed: {}", e));
        register_slot(0, &ref_bytes)
            .unwrap_or_else(|e| panic!("CUDA: failed to upload reference set: {}", e));
        if let Some(eb) = excl_bytes.as_ref() {
            register_slot(1, eb)
                .unwrap_or_else(|e| panic!("CUDA: failed to upload exclusivity set: {}", e));
        }
    }

    let lengths = build_length_list(
        params.min_oligo_length,
        params.max_oligo_length,
        params.length_skip,
    );
    let total_lengths = lengths.len() as u32;

    for (length_idx, oligo_length) in lengths.into_iter().enumerate() {
        let ref_bytes = Arc::clone(&ref_bytes);
        let excl_bytes = excl_bytes.clone();
        let excl_names = excl_names.clone();
        let length_result = pool.install(|| {
            analyze_length(
                template,
                &ref_bytes,
                excl_bytes.as_ref().map(|v| v.as_slice()),
                excl_names.as_ref().map(|v| v.as_slice()),
                params,
                oligo_length,
                length_idx as u32,
                total_lengths,
                &progress_tx,
                pause.as_ref(),
            )
        });

        results.results_by_length.insert(oligo_length, length_result);
    }

    results
}

fn analyze_length(
    template: &TemplateData,
    ref_bytes: &[Vec<u8>],
    excl_bytes: Option<&[Vec<u8>]>,
    excl_names: Option<&[String]>,
    params: &AnalysisParams,
    oligo_length: u32,
    length_idx: u32,
    total_lengths: u32,
    progress_tx: &Option<Sender<ProgressUpdate>>,
    pause: Option<&PauseFlag>,
) -> LengthResult {
    let length = oligo_length as usize;
    let resolution = params.resolution as usize;
    let template_len = template.sequence.len();

    let max_start = if template_len >= length { template_len - length } else { 0 };
    let positions: Vec<usize> = (0..=max_start).step_by(resolution).collect();
    let total_positions = positions.len();

    let completed_count = Arc::new(AtomicUsize::new(0));
    let template_bytes = template.sequence.as_bytes();

    let max_ref_len = ref_bytes.iter().map(|r| r.len()).max().unwrap_or(0);
    let max_excl_len = excl_bytes
        .map(|eb| eb.iter().map(|r| r.len()).max().unwrap_or(0))
        .unwrap_or(0);
    let max_seq_len = max_ref_len.max(max_excl_len);
    let pw_params = params.pairwise;
    let aligner_kind = params.aligner;
    let simple_params = params.simple.unwrap_or_default();

    let mut position_results: Vec<PositionResult> = positions
        .par_iter()
        .map_init(
            move || match aligner_kind {
                AlignerKind::Pairwise => {
                    WorkerAligner::Pairwise(create_aligner(length, max_seq_len, &pw_params))
                }
                AlignerKind::Simple => {
                    WorkerAligner::Simple(create_simple_aligner(&simple_params))
                }
                AlignerKind::SimpleSimd => {
                    WorkerAligner::SimpleSimd(create_simd_aligner(&simple_params))
                }
                #[cfg(feature = "cuda")]
                AlignerKind::SimpleCuda => {
                    WorkerAligner::SimpleCuda(create_cuda_aligner(&simple_params))
                }
                #[cfg(not(feature = "cuda"))]
                AlignerKind::SimpleCuda => {
                    panic!(
                        "aligner = simple_cuda was selected but this build was \
                         compiled without the `cuda` feature — rebuild with \
                         `cargo build --features cuda`"
                    );
                }
            },
            |aligner, &position| {
                if let Some(p) = pause {
                    p.wait_if_paused();
                }
                let analysis = analyze_window(
                    template_bytes,
                    ref_bytes,
                    params,
                    position,
                    length,
                    aligner,
                );

                let exclusivity = excl_bytes.map(|eb| {
                    analyze_exclusivity(
                        template_bytes,
                        eb,
                        excl_names.unwrap(),
                        params,
                        position,
                        length,
                        aligner,
                    )
                });

                let completed = completed_count.fetch_add(1, Ordering::Relaxed) + 1;
                if let Some(tx) = progress_tx {
                    if completed % 10 == 0 || completed == total_positions {
                        let _ = tx.send(ProgressUpdate {
                            current_length: oligo_length,
                            current_position: position,
                            total_positions,
                            lengths_completed: length_idx,
                            total_lengths,
                            message: format!(
                                "Length {}/{}: Position {}/{}",
                                length_idx + 1,
                                total_lengths,
                                completed,
                                total_positions
                            ),
                        });
                    }
                }

                PositionResult {
                    position,
                    variants_needed: analysis.variants_for_threshold,
                    analysis,
                    exclusivity,
                }
            },
        )
        .collect();

    position_results.sort_by_key(|r| r.position);

    LengthResult {
        oligo_length,
        positions: position_results,
    }
}

fn analyze_window(
    template_bytes: &[u8],
    ref_bytes: &[Vec<u8>],
    params: &AnalysisParams,
    position: usize,
    length: usize,
    aligner: &mut WorkerAligner,
) -> WindowAnalysisResult {
    let oligo = &template_bytes[position..position + length];
    let total_refs = ref_bytes.len();

    let (matched_sequences, no_match_count) = match aligner {
        WorkerAligner::Pairwise(a) => {
            collect_matches_with_aligner(a, oligo, ref_bytes, &params.pairwise)
        }
        WorkerAligner::Simple(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_matches_with_simple_aligner(a, oligo, ref_bytes, &sp)
        }
        WorkerAligner::SimpleSimd(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_matches_with_simd_aligner(a, oligo, ref_bytes, &sp)
        }
        #[cfg(feature = "cuda")]
        WorkerAligner::SimpleCuda(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_matches_with_cuda_aligner(a, oligo, ref_bytes, &sp)
        }
    };

    if matched_sequences.is_empty() {
        return WindowAnalysisResult {
            total_sequences: total_refs,
            sequences_analyzed: 0,
            no_match_count,
            skipped: true,
            skip_reason: Some("No valid matches found in any reference sequence".to_string()),
            ..Default::default()
        };
    }

    let seq_refs: Vec<&str> = matched_sequences.iter().map(|s| s.as_str()).collect();

    let mut result = analyze_sequences(
        &seq_refs,
        &params.method,
        params.exclude_n,
        params.coverage_threshold,
        params.max_seeds as usize,
    );

    result.total_sequences = total_refs;
    result.sequences_analyzed = matched_sequences.len();
    result.no_match_count = no_match_count;

    if total_refs > matched_sequences.len() {
        let total_f = total_refs as f64;
        for variant in &mut result.variants {
            variant.percentage = (variant.count as f64 / total_f) * 100.0;
        }
    }

    let var_limit_applied = if let Some(limit) = params.var_limit {
        let limit = limit as usize;
        if result.variants.len() > limit {
            let dropped_count: usize = result.variants[limit..].iter().map(|v| v.count).sum();
            result.variants.truncate(limit);
            result.no_match_count += dropped_count;
            result.sequences_analyzed = result.sequences_analyzed.saturating_sub(dropped_count);
            true
        } else {
            false
        }
    } else {
        false
    };

    if total_refs > matched_sequences.len() || var_limit_applied {
        let mut cumulative = 0.0;
        let mut new_variants_needed = result.variants.len();
        let mut new_coverage = 0.0;
        for (i, variant) in result.variants.iter().enumerate() {
            cumulative += variant.percentage;
            if cumulative >= params.coverage_threshold {
                new_variants_needed = i + 1;
                new_coverage = cumulative;
                break;
            }
        }
        if cumulative < params.coverage_threshold {
            new_coverage = cumulative;
        }
        result.variants_for_threshold = new_variants_needed;
        result.coverage_at_threshold = new_coverage;
    }

    result
}

fn analyze_exclusivity(
    template_bytes: &[u8],
    excl_bytes: &[Vec<u8>],
    excl_names: &[String],
    params: &AnalysisParams,
    position: usize,
    length: usize,
    aligner: &mut WorkerAligner,
) -> ExclusivityResult {
    let oligo = &template_bytes[position..position + length];
    let mismatch_counts = match aligner {
        WorkerAligner::Pairwise(a) => {
            collect_mismatch_counts_with_aligner(a, oligo, excl_bytes, &params.pairwise)
        }
        WorkerAligner::Simple(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_mismatch_counts_with_simple_aligner(a, oligo, excl_bytes, &sp)
        }
        WorkerAligner::SimpleSimd(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_mismatch_counts_with_simd_aligner(a, oligo, excl_bytes, &sp)
        }
        #[cfg(feature = "cuda")]
        WorkerAligner::SimpleCuda(a) => {
            let sp = params.simple.unwrap_or_default();
            collect_mismatch_counts_with_cuda_aligner(a, oligo, excl_bytes, &sp)
        }
    };

    let mut buckets: std::collections::HashMap<u32, (usize, String)> =
        std::collections::HashMap::new();
    let mut no_match_count = 0usize;
    let mut no_match_example = String::new();
    let mut min_mismatches: Option<u32> = None;

    for (i, mm) in mismatch_counts.iter().enumerate() {
        match mm {
            Some(m) => {
                let entry = buckets.entry(*m).or_insert_with(|| (0, excl_names[i].clone()));
                entry.0 += 1;
                match min_mismatches {
                    None => min_mismatches = Some(*m),
                    Some(current) if *m < current => min_mismatches = Some(*m),
                    _ => {}
                }
            }
            None => {
                if no_match_count == 0 {
                    no_match_example = excl_names[i].clone();
                }
                no_match_count += 1;
            }
        }
    }

    let mut mismatch_histogram: Vec<MismatchBucket> = buckets
        .into_iter()
        .map(|(mismatches, (count, example_name))| MismatchBucket {
            mismatches,
            count,
            example_name,
        })
        .collect();
    mismatch_histogram.sort_by_key(|b| b.mismatches);

    if no_match_count > 0 {
        mismatch_histogram.push(MismatchBucket {
            mismatches: u32::MAX,
            count: no_match_count,
            example_name: no_match_example,
        });
    }

    ExclusivityResult {
        total_sequences: excl_bytes.len(),
        no_match_count,
        mismatch_histogram,
        min_mismatches,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::types::AnalysisMethod;

    #[test]
    fn test_build_length_list() {
        assert_eq!(build_length_list(20, 25, 0), vec![20, 21, 22, 23, 24, 25]);
        assert_eq!(build_length_list(20, 25, 1), vec![20, 22, 24]);
        assert_eq!(build_length_list(20, 25, 2), vec![20, 23]);
        assert_eq!(build_length_list(18, 18, 0), vec![18]);
        assert_eq!(build_length_list(18, 18, 5), vec![18]);
        assert!(build_length_list(25, 20, 0).is_empty());
    }

    #[test]
    fn test_screening_example() {
        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
        };

        let references = ReferenceData {
            names: vec![
                "Ref1".to_string(),
                "Ref2".to_string(),
                "Ref3".to_string(),
                "Ref4".to_string(),
            ],
            sequences: vec![
                "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_string(),
                "GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
            ],
        };

        let params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            min_oligo_length: 10,
            max_oligo_length: 10,
            resolution: 1,
            coverage_threshold: 95.0,
            ..Default::default()
        };

        let results = run_screening(&template, &references, &params, None, None, None);
        assert!(results.results_by_length.contains_key(&10));

        let length_result = results.results_by_length.get(&10).unwrap();
        let first_pos = &length_result.positions[0];
        assert!(!first_pos.analysis.skipped);
        assert!(!first_pos.analysis.variants.is_empty());
        assert!(first_pos.exclusivity.is_none());
    }

    #[test]
    fn test_screening_with_exclusivity() {
        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
        };

        let references = ReferenceData {
            names: vec!["Ref1".to_string()],
            sequences: vec!["TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string()],
        };

        let exclusivity = ReferenceData {
            names: vec!["Excl1".to_string(), "Excl2".to_string()],
            sequences: vec![
                "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string(),
            ],
        };

        let params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            min_oligo_length: 10,
            max_oligo_length: 10,
            resolution: 1,
            coverage_threshold: 95.0,
            ..Default::default()
        };

        let results = run_screening(&template, &references, &params, Some(&exclusivity), None, None);
        let length_result = results.results_by_length.get(&10).unwrap();
        let first_pos = &length_result.positions[0];

        assert!(first_pos.exclusivity.is_some());
        let excl = first_pos.exclusivity.as_ref().unwrap();
        assert_eq!(excl.total_sequences, 2);
        assert!(results.differential_enabled);
        assert_eq!(results.exclusivity_sequence_count, Some(2));
    }

    #[test]
    fn var_limit_folds_dropped_variants_into_no_match() {
        // Use the simple (bitap) aligner with max_mismatches=2 so each
        // single-base-substituted reference returns a full-length matched
        // fragment. Without a limit this position would have 4 variants;
        // with var_limit = 2 the bottom two are dropped and folded into
        // no_match.
        use crate::analysis::types::{AlignerKind, SimpleParams};

        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "ACACACACACACACACACAC".to_string(),
        };
        let references = ReferenceData {
            names: vec!["R1".into(), "R2".into(), "R3".into(), "R4".into()],
            sequences: vec![
                "ACACACACACACACACACAC".into(),
                "ACACACACACACACACACAA".into(),
                "ACACACACACACACACACAG".into(),
                "ACACACACACACACACACAT".into(),
            ],
        };

        let params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            aligner: AlignerKind::Simple,
            simple: Some(SimpleParams { max_mismatches: 2 }),
            min_oligo_length: 20,
            max_oligo_length: 20,
            resolution: 1,
            coverage_threshold: 95.0,
            var_limit: Some(2),
            ..Default::default()
        };

        let results = run_screening(&template, &references, &params, None, None, None);
        let pos = &results.results_by_length.get(&20).unwrap().positions[0];
        let a = &pos.analysis;

        assert_eq!(a.variants.len(), 2);
        assert_eq!(a.total_sequences, 4);
        assert_eq!(a.sequences_analyzed, 2);
        assert_eq!(a.no_match_count, 2);
        let kept_count: usize = a.variants.iter().map(|v| v.count).sum();
        assert_eq!(kept_count + a.no_match_count, a.total_sequences);
    }

    /// SimpleSimd is required to produce bit-identical results to Simple on
    /// the same inputs (same recurrence, same per-reference state). This is
    /// the regression test that pins that invariant; if you see it fail,
    /// the SIMD kernel diverged from the scalar bitap.
    #[test]
    fn simple_simd_matches_simple_on_screening_results() {
        use crate::analysis::types::{AlignerKind, SimpleParams};

        #[cfg(target_arch = "x86_64")]
        {
            if !std::is_x86_feature_detected!("avx2") {
                eprintln!("skipping simple_simd equivalence test: AVX2 not available");
                return;
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            eprintln!("skipping simple_simd equivalence test: not x86_64");
            return;
        }

        // Template + references chosen to exercise: full SIMD batches,
        // a non-LANES-divisible tail, forward and reverse-complement matches,
        // mismatches up to k, and references shorter than the longest in
        // their batch.
        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "TATGGTACGTCATGTTCTAGAAATGGGCTGTACGTACGTACGTACGTACGTACGT".to_string(),
        };
        let references = ReferenceData {
            names: (0..11).map(|i| format!("R{}", i)).collect(),
            sequences: vec![
                "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),     // exact fwd
                "AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),   // exact fwd, offset
                "TATGGTTCGTCATGTTCTAGAAATGGGCTGT".to_string(),     // 1mm fwd
                "GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),    // exact fwd, offset
                // RC of "TATGGTACGT" is "ACGTACCATA"
                "NNNNACGTACCATANNNNNNNNNNNNNNNNNN".to_string(),    // RC match
                "GGGGGGGGGGGG".to_string(),                        // shorter, no match
                "TATGGTACGTCATG".to_string(),                      // exact fwd, short ref
                "NNNNNNNN".to_string(),                            // too short / no match
                "TATGGTACGGCATGTTCT".to_string(),                  // 1mm
                "TATGGTACGT".to_string(),                          // exact, length == oligo
                "TATGGTACATCATG".to_string(),                      // 1mm fwd
            ],
        };
        let exclusivity = ReferenceData {
            names: (0..6).map(|i| format!("E{}", i)).collect(),
            sequences: vec![
                "TATGGTACGTCATG".to_string(),     // exact
                "AAAAAAAAAAAAAA".to_string(),     // none
                "TATGGTTCGTCATG".to_string(),     // 1 mm
                "NNNNACGTACCATANNNN".to_string(), // RC
                "GG".to_string(),                 // too short
                "TATGGTACATCATG".to_string(),     // 1 mm
            ],
        };

        let mut params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            aligner: AlignerKind::Simple,
            simple: Some(SimpleParams { max_mismatches: 2 }),
            min_oligo_length: 10,
            max_oligo_length: 14,
            resolution: 1,
            coverage_threshold: 95.0,
            ..Default::default()
        };

        let simple_results = run_screening(
            &template,
            &references,
            &params,
            Some(&exclusivity),
            None,
            None,
        );

        params.aligner = AlignerKind::SimpleSimd;
        let simd_results = run_screening(
            &template,
            &references,
            &params,
            Some(&exclusivity),
            None,
            None,
        );

        // Normalize before comparing:
        //   1. Drop the `aligner` field (the only legitimate difference).
        //   2. Sort every `variants` array by (count desc, sequence asc).
        //      `find_variants_no_ambiguities` builds variants from a HashMap
        //      and stable-sorts by count only; tied counts keep nondeterministic
        //      HashMap iteration order, which would otherwise cause two runs
        //      of the *same* aligner to diverge in serialization. The SIMD
        //      backend is exact-equivalent to scalar at the matched-fragment
        //      level — this normalization isolates that invariant.
        fn normalize(v: &mut serde_json::Value) {
            match v {
                serde_json::Value::Object(map) => {
                    map.remove("aligner");
                    for (k, child) in map.iter_mut() {
                        if k == "variants" {
                            if let Some(arr) = child.as_array_mut() {
                                arr.sort_by(|a, b| {
                                    let ca = a.get("count").and_then(|x| x.as_u64()).unwrap_or(0);
                                    let cb = b.get("count").and_then(|x| x.as_u64()).unwrap_or(0);
                                    cb.cmp(&ca).then_with(|| {
                                        let sa = a
                                            .get("sequence")
                                            .and_then(|x| x.as_str())
                                            .unwrap_or("");
                                        let sb = b
                                            .get("sequence")
                                            .and_then(|x| x.as_str())
                                            .unwrap_or("");
                                        sa.cmp(sb)
                                    })
                                });
                            }
                        }
                        normalize(child);
                    }
                }
                serde_json::Value::Array(arr) => {
                    for child in arr.iter_mut() {
                        normalize(child);
                    }
                }
                _ => {}
            }
        }
        let mut lhs = serde_json::to_value(&simple_results).unwrap();
        let mut rhs = serde_json::to_value(&simd_results).unwrap();
        normalize(&mut lhs);
        normalize(&mut rhs);
        if lhs != rhs {
            let lhs_s = serde_json::to_string_pretty(&lhs).unwrap();
            let rhs_s = serde_json::to_string_pretty(&rhs).unwrap();
            let lhs_lines: Vec<&str> = lhs_s.lines().collect();
            let rhs_lines: Vec<&str> = rhs_s.lines().collect();
            for (i, (l, r)) in lhs_lines.iter().zip(rhs_lines.iter()).enumerate() {
                if l != r {
                    let lo = i.saturating_sub(3);
                    let hi = (i + 4).min(lhs_lines.len()).min(rhs_lines.len());
                    eprintln!("first diff context (lines {}..{}):", lo + 1, hi);
                    for j in lo..hi {
                        let mark = if j == i { ">>" } else { "  " };
                        eprintln!("{} L{:>3}: simple      | {}", mark, j + 1, lhs_lines[j]);
                        eprintln!("{} L{:>3}: simple_simd | {}", mark, j + 1, rhs_lines[j]);
                    }
                    break;
                }
            }
            panic!("simple_simd diverged from simple");
        }
    }

    /// Bit-identical equivalence between the scalar bitap (`simple`) and the
    /// CUDA bitap (`simple_cuda`) on a fixture that exercises full SIMD/CUDA
    /// batches, ragged-length tails, forward + reverse-complement matches,
    /// and exclusivity scanning. Mirrors `simple_simd_matches_simple_on_screening_results`.
    ///
    /// Uses `--method none` only — and that is *sufficient* to cover every
    /// method. The aligner backend is responsible solely for the alignment
    /// stage (`collect_matches_with_*` / `collect_mismatch_counts_with_*`);
    /// `none` compares that stage's output directly. The downstream variant
    /// analysis (`analyze_sequences`) is backend-agnostic — it consumes the
    /// matched fragments and never sees the aligner — so once the fragments
    /// match, `fixed` and `incremental` follow by construction.
    ///
    /// Skips if no CUDA-capable GPU is available (build env can compile the
    /// feature without a runtime GPU).
    #[cfg(feature = "cuda")]
    #[test]
    fn simple_cuda_matches_simple_on_screening_results() {
        use crate::analysis::cuda_align::ensure_initialized;
        use crate::analysis::types::{AlignerKind, SimpleParams};

        if ensure_initialized().is_err() {
            eprintln!("skipping simple_cuda equivalence test: no usable CUDA device");
            return;
        }

        let template = TemplateData {
            name: "Template".to_string(),
            sequence: "TATGGTACGTCATGTTCTAGAAATGGGCTGTACGTACGTACGTACGTACGTACGT".to_string(),
        };
        let references = ReferenceData {
            names: (0..11).map(|i| format!("R{}", i)).collect(),
            sequences: vec![
                "TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "TATGGTTCGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_string(),
                "NNNNACGTACCATANNNNNNNNNNNNNNNNNN".to_string(),
                "GGGGGGGGGGGG".to_string(),
                "TATGGTACGTCATG".to_string(),
                "NNNNNNNN".to_string(),
                "TATGGTACGGCATGTTCT".to_string(),
                "TATGGTACGT".to_string(),
                "TATGGTACATCATG".to_string(),
            ],
        };
        let exclusivity = ReferenceData {
            names: (0..6).map(|i| format!("E{}", i)).collect(),
            sequences: vec![
                "TATGGTACGTCATG".to_string(),
                "AAAAAAAAAAAAAA".to_string(),
                "TATGGTTCGTCATG".to_string(),
                "NNNNACGTACCATANNNN".to_string(),
                "GG".to_string(),
                "TATGGTACATCATG".to_string(),
            ],
        };

        let mut params = AnalysisParams {
            method: AnalysisMethod::NoAmbiguities,
            aligner: AlignerKind::Simple,
            simple: Some(SimpleParams { max_mismatches: 2 }),
            min_oligo_length: 10,
            max_oligo_length: 14,
            resolution: 1,
            coverage_threshold: 95.0,
            ..Default::default()
        };

        let simple_results = run_screening(
            &template,
            &references,
            &params,
            Some(&exclusivity),
            None,
            None,
        );

        params.aligner = AlignerKind::SimpleCuda;
        let cuda_results = run_screening(
            &template,
            &references,
            &params,
            Some(&exclusivity),
            None,
            None,
        );

        // Same normalization as the SIMD equivalence test.
        fn normalize(v: &mut serde_json::Value) {
            match v {
                serde_json::Value::Object(map) => {
                    map.remove("aligner");
                    for (k, child) in map.iter_mut() {
                        if k == "variants" {
                            if let Some(arr) = child.as_array_mut() {
                                arr.sort_by(|a, b| {
                                    let ca = a.get("count").and_then(|x| x.as_u64()).unwrap_or(0);
                                    let cb = b.get("count").and_then(|x| x.as_u64()).unwrap_or(0);
                                    cb.cmp(&ca).then_with(|| {
                                        let sa = a
                                            .get("sequence")
                                            .and_then(|x| x.as_str())
                                            .unwrap_or("");
                                        let sb = b
                                            .get("sequence")
                                            .and_then(|x| x.as_str())
                                            .unwrap_or("");
                                        sa.cmp(sb)
                                    })
                                });
                            }
                        }
                        normalize(child);
                    }
                }
                serde_json::Value::Array(arr) => {
                    for child in arr.iter_mut() {
                        normalize(child);
                    }
                }
                _ => {}
            }
        }
        let mut lhs = serde_json::to_value(&simple_results).unwrap();
        let mut rhs = serde_json::to_value(&cuda_results).unwrap();
        normalize(&mut lhs);
        normalize(&mut rhs);
        if lhs != rhs {
            let lhs_s = serde_json::to_string_pretty(&lhs).unwrap();
            let rhs_s = serde_json::to_string_pretty(&rhs).unwrap();
            let lhs_lines: Vec<&str> = lhs_s.lines().collect();
            let rhs_lines: Vec<&str> = rhs_s.lines().collect();
            for (i, (l, r)) in lhs_lines.iter().zip(rhs_lines.iter()).enumerate() {
                if l != r {
                    let lo = i.saturating_sub(3);
                    let hi = (i + 4).min(lhs_lines.len()).min(rhs_lines.len());
                    eprintln!("first diff context (lines {}..{}):", lo + 1, hi);
                    for j in lo..hi {
                        let mark = if j == i { ">>" } else { "  " };
                        eprintln!("{} L{:>3}: simple      | {}", mark, j + 1, lhs_lines[j]);
                        eprintln!("{} L{:>3}: simple_cuda | {}", mark, j + 1, rhs_lines[j]);
                    }
                    break;
                }
            }
            panic!("simple_cuda diverged from simple");
        }
    }
}
