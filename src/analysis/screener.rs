//! Screening logic using pairwise alignment.
//!
//! Iterates through the template sequence with different oligo lengths,
//! using pairwise alignment to find best matches in each reference sequence.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::Sender;

use rayon::prelude::*;

use super::analyzer::analyze_sequences;
use super::fasta::{ReferenceData, TemplateData};
use super::pairwise::{
    DnaAligner, collect_matches_with_aligner, collect_mismatch_counts_with_aligner, create_aligner,
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
    );

    result.total_sequences = total_refs;
    result.sequences_analyzed = matched_sequences.len();
    result.no_match_count = no_match_count;

    if total_refs > matched_sequences.len() {
        let total_f = total_refs as f64;
        for variant in &mut result.variants {
            variant.percentage = (variant.count as f64 / total_f) * 100.0;
        }
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
}
