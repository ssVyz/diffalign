//! AVX2-vectorized bitap aligner backend.
//!
//! Mirrors the API surface of `simple_align.rs` so it plugs into the
//! `WorkerAligner` dispatch in `screener.rs` without further changes:
//!
//! * `SimdAligner` — per-worker reusable state.
//! * `create_simd_aligner` — constructed once via Rayon's `map_init`.
//! * `collect_matches_with_simd_aligner` / `collect_mismatch_counts_with_simd_aligner`
//!
//! Internally, references are processed in chunks of `LANES_AVX2` (4) using
//! the AVX2 kernel; any trailing references that don't fill a full chunk are
//! processed with the scalar bitap. Forward and reverse-complement strands are
//! both scanned, with the same tiebreaker as `simple_align::best_hit`
//! (lowest mismatches → earliest start → forward preferred over reverse).
//!
//! The output is bit-identical to `simple_align` on the same inputs: bitap is
//! exact integer arithmetic and each lane's state evolves independently of the
//! others. The equivalence test in `screener.rs` enforces this invariant.

use super::simplescreen::bitap::BitapState;
use super::simplescreen::pattern::{MAX_PATTERN_LEN, PreparedPattern};
use super::simplescreen::screener::{Hit, Orientation, screen};
use super::types::SimpleParams;

#[cfg(target_arch = "x86_64")]
use super::simplescreen::bitap_simd::{LANES_AVX2, LaneBest, scan_avx2};

#[cfg(target_arch = "x86_64")]
const LANES: usize = LANES_AVX2;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 1;

/// Per-worker scratch state for the SIMD aligner.
///
/// `Send`, not `Sync`. The internal scalar `BitapState` + hit buffer are used
/// for the tail of each reference chunk (whatever doesn't fill a full SIMD
/// batch) and as a fallback path on non-x86 builds.
pub struct SimdAligner {
    scalar_state: BitapState,
    scalar_hits: Vec<Hit>,
}

pub fn create_simd_aligner(params: &SimpleParams) -> SimdAligner {
    SimdAligner {
        scalar_state: BitapState::new(params.max_mismatches),
        scalar_hits: Vec::with_capacity(16),
    }
}

/// Largest oligo length the SIMD backend can handle. Same single-u64 cap as
/// the scalar bitap.
pub const SIMD_MAX_OLIGO_LEN: usize = MAX_PATTERN_LEN;

/// Same shape as `simple_align::SimpleMatch`. Kept private to avoid leaking a
/// duplicate name; the public API returns matched fragments as `String` just
/// like `simple_align`.
#[derive(Debug, Clone)]
struct SimdMatch {
    matched_sequence: String,
    mismatches: u32,
    orientation: Orientation,
}

fn complement_byte(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' | b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

fn reverse_complement(bytes: &[u8]) -> String {
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        out.push(complement_byte(b));
    }
    String::from_utf8(out).unwrap_or_default()
}

/// Run the scalar bitap on one reference. Used for the tail of a reference
/// slice that doesn't fill a full SIMD batch, and as the implementation on
/// non-x86 builds (never reached in practice — `into_analysis_params` rejects
/// `simple_simd` on non-x86 targets).
fn scan_reference_scalar(
    aligner: &mut SimdAligner,
    pattern: &PreparedPattern,
    reference: &[u8],
) -> Option<SimdMatch> {
    aligner.scalar_hits.clear();
    screen(&mut aligner.scalar_state, pattern, reference, &mut aligner.scalar_hits);

    let hit = best_scalar_hit(&aligner.scalar_hits)?;
    let len = pattern.len as usize;
    let start_0 = (hit.start as usize) - 1;
    let end_0 = hit.end as usize;
    let frag = &reference[start_0..end_0];
    debug_assert_eq!(frag.len(), len);

    let matched_sequence = match hit.orientation {
        Orientation::Forward => String::from_utf8_lossy(frag).into_owned(),
        Orientation::Reverse => reverse_complement(frag),
    };

    Some(SimdMatch {
        matched_sequence,
        mismatches: hit.mismatches,
        orientation: hit.orientation,
    })
}

/// Tiebreaker identical to `simple_align::best_hit`.
fn best_scalar_hit(hits: &[Hit]) -> Option<Hit> {
    hits.iter()
        .min_by(|a, b| {
            a.mismatches
                .cmp(&b.mismatches)
                .then_with(|| a.start.cmp(&b.start))
                .then_with(|| match (a.orientation, b.orientation) {
                    (Orientation::Forward, Orientation::Reverse) => std::cmp::Ordering::Less,
                    (Orientation::Reverse, Orientation::Forward) => std::cmp::Ordering::Greater,
                    _ => std::cmp::Ordering::Equal,
                })
        })
        .copied()
}

#[cfg(target_arch = "x86_64")]
fn scan_batch_avx2(
    pattern: &PreparedPattern,
    batch: &[&[u8]; LANES_AVX2],
) -> [(Option<u32>, Option<u32>, Orientation); LANES_AVX2] {
    // (best_mm, best_end_pos, orientation) per lane; None if no hit.
    let k = pattern.max_mismatches;
    let end_bit = pattern.end_bit;
    let fwd = unsafe { scan_avx2(&pattern.fwd_masks, batch, k, end_bit) };
    let rev = unsafe { scan_avx2(&pattern.rc_masks, batch, k, end_bit) };

    let mut out: [(Option<u32>, Option<u32>, Orientation); LANES_AVX2] = [
        (None, None, Orientation::Forward),
        (None, None, Orientation::Forward),
        (None, None, Orientation::Forward),
        (None, None, Orientation::Forward),
    ];
    for lane in 0..LANES_AVX2 {
        out[lane] = combine_fwd_rev(&fwd[lane], &rev[lane]);
    }
    out
}

#[cfg(target_arch = "x86_64")]
fn combine_fwd_rev(
    fwd: &LaneBest,
    rev: &LaneBest,
) -> (Option<u32>, Option<u32>, Orientation) {
    match (fwd.found, rev.found) {
        (false, false) => (None, None, Orientation::Forward),
        (true, false) => (Some(fwd.mismatches), Some(fwd.end_pos), Orientation::Forward),
        (false, true) => (Some(rev.mismatches), Some(rev.end_pos), Orientation::Reverse),
        (true, true) => {
            // Tiebreaker: lower mismatches; on equal mm, lower end_pos
            // (== lower start, since pattern length is fixed); on equal mm and
            // end, forward wins.
            if fwd.mismatches < rev.mismatches
                || (fwd.mismatches == rev.mismatches && fwd.end_pos <= rev.end_pos)
            {
                (Some(fwd.mismatches), Some(fwd.end_pos), Orientation::Forward)
            } else {
                (Some(rev.mismatches), Some(rev.end_pos), Orientation::Reverse)
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
fn extract_simd_match(
    pattern: &PreparedPattern,
    reference: &[u8],
    lane_result: (Option<u32>, Option<u32>, Orientation),
) -> Option<SimdMatch> {
    let (mm, end_pos, orientation) = lane_result;
    let mm = mm?;
    let end_pos = end_pos?;
    let len = pattern.len as usize;
    let end_0 = end_pos as usize;
    let start_0 = end_0 + 1 - len;
    let frag = &reference[start_0..end_0 + 1];
    debug_assert_eq!(frag.len(), len);
    let matched_sequence = match orientation {
        Orientation::Forward => String::from_utf8_lossy(frag).into_owned(),
        Orientation::Reverse => reverse_complement(frag),
    };
    Some(SimdMatch {
        matched_sequence,
        mismatches: mm,
        orientation,
    })
}

pub fn collect_matches_with_simd_aligner(
    aligner: &mut SimdAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> (Vec<String>, usize) {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return (Vec::new(), references.len()),
    };
    if !pattern.valid {
        return (Vec::new(), references.len());
    }

    let mut matched = Vec::with_capacity(references.len());
    let mut no_match_count = 0usize;

    // SIMD batches (x86_64). On non-x86 this loop body is a no-op; everything
    // falls through to the scalar tail loop below.
    #[cfg(target_arch = "x86_64")]
    {
        let mut idx = 0;
        while idx + LANES <= references.len() {
            let batch: [&[u8]; LANES_AVX2] = [
                references[idx].as_slice(),
                references[idx + 1].as_slice(),
                references[idx + 2].as_slice(),
                references[idx + 3].as_slice(),
            ];
            let results = scan_batch_avx2(&pattern, &batch);
            for lane in 0..LANES_AVX2 {
                match extract_simd_match(&pattern, references[idx + lane].as_slice(), results[lane])
                {
                    Some(m) => matched.push(m.matched_sequence),
                    None => no_match_count += 1,
                }
            }
            idx += LANES;
        }
        // Tail: < LANES references left.
        for r in &references[idx..] {
            match scan_reference_scalar(aligner, &pattern, r) {
                Some(m) => matched.push(m.matched_sequence),
                None => no_match_count += 1,
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        for r in references {
            match scan_reference_scalar(aligner, &pattern, r) {
                Some(m) => matched.push(m.matched_sequence),
                None => no_match_count += 1,
            }
        }
    }

    (matched, no_match_count)
}

pub fn collect_mismatch_counts_with_simd_aligner(
    aligner: &mut SimdAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> Vec<Option<u32>> {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return references.iter().map(|_| None).collect(),
    };
    if !pattern.valid {
        return references.iter().map(|_| None).collect();
    }

    let mut out: Vec<Option<u32>> = Vec::with_capacity(references.len());

    #[cfg(target_arch = "x86_64")]
    {
        let mut idx = 0;
        while idx + LANES <= references.len() {
            let batch: [&[u8]; LANES_AVX2] = [
                references[idx].as_slice(),
                references[idx + 1].as_slice(),
                references[idx + 2].as_slice(),
                references[idx + 3].as_slice(),
            ];
            let results = scan_batch_avx2(&pattern, &batch);
            for lane in 0..LANES_AVX2 {
                out.push(results[lane].0);
            }
            idx += LANES;
        }
        for r in &references[idx..] {
            out.push(scan_reference_scalar(aligner, &pattern, r).map(|m| m.mismatches));
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        for r in references {
            out.push(scan_reference_scalar(aligner, &pattern, r).map(|m| m.mismatches));
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> SimpleParams {
        SimpleParams { max_mismatches: 2 }
    }

    fn skip_if_no_simd() -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                return false;
            }
        }
        eprintln!("skipping SIMD-aligner test: AVX2 not available");
        true
    }

    #[test]
    fn forward_exact_match_via_simd() {
        if skip_if_no_simd() {
            return;
        }
        let mut aligner = create_simd_aligner(&p());
        let oligo = b"TATGGTACGT";
        let refs: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"NNNNNNNNNNNNNNNNNNNNN".to_vec(),
            b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_vec(),
        ];
        let (matched, nm) = collect_matches_with_simd_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 1);
        // Order-preserving: matched[0] from refs[0], etc.
        assert_eq!(matched, vec![
            "TATGGTACGT".to_string(),
            "TATGGTACGT".to_string(),
            "TATGGTTCGT".to_string(),
        ]);
    }

    #[test]
    fn reverse_complement_returned_in_oligo_orientation() {
        if skip_if_no_simd() {
            return;
        }
        let mut aligner = create_simd_aligner(&p());
        let oligo = b"TATGGTACGT";
        // RC of oligo = ACGTACCATA.
        let refs: Vec<Vec<u8>> = vec![b"NNNNACGTACCATANNNN".to_vec()];
        let (matched, nm) = collect_matches_with_simd_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 0);
        assert_eq!(matched, vec!["TATGGTACGT".to_string()]);
    }

    #[test]
    fn tail_chunk_falls_back_to_scalar() {
        if skip_if_no_simd() {
            return;
        }
        // 5 refs → one SIMD batch of 4 + one scalar tail.
        let mut aligner = create_simd_aligner(&p());
        let oligo = b"ACGT";
        let refs: Vec<Vec<u8>> = vec![
            b"NNACGT".to_vec(),
            b"NNACGT".to_vec(),
            b"NNACGT".to_vec(),
            b"NNACGT".to_vec(),
            b"NNACGT".to_vec(),
        ];
        let (matched, nm) = collect_matches_with_simd_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 0);
        assert_eq!(matched.len(), 5);
        for m in &matched {
            assert_eq!(m, "ACGT");
        }
    }

    #[test]
    fn mismatch_counts_match_inputs_per_reference() {
        if skip_if_no_simd() {
            return;
        }
        let mut aligner = create_simd_aligner(&p());
        let oligo = b"TATGGTACGT";
        let refs: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATG".to_vec(), // exact
            b"GGGGGGGGGGGGGG".to_vec(), // none
            b"TATGGTTCGTCATG".to_vec(), // 1 mm
            b"TATGGGTCGTCATG".to_vec(), // 2 mm
            b"TATGGTACGTNNNN".to_vec(), // exact (tail)
        ];
        let counts = collect_mismatch_counts_with_simd_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(counts.len(), 5);
        assert_eq!(counts[0], Some(0));
        assert_eq!(counts[1], None);
        assert_eq!(counts[2], Some(1));
        assert_eq!(counts[3], Some(2));
        assert_eq!(counts[4], Some(0));
    }
}
