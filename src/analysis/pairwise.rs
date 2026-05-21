//! Pairwise alignment for matching oligos against reference sequences.
//!
//! Forward + reverse-complement scanning is enabled. When the best hit for a
//! reference is on the reverse strand, the matched fragment is
//! reverse-complemented before being returned, so downstream variant analysis
//! sees all fragments in the same orientation as the oligo. Mirrors the
//! strand-handling contract of the bitap aligners in `simple_align.rs`.

use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::{Aligner, MatchFunc, MatchParams};

use super::simplescreen::screener::Orientation;
use super::types::{AnchorHit, PairwiseParams};

pub type DnaAligner = Aligner<MatchParams>;

pub fn create_aligner(oligo_len: usize, max_ref_len: usize, params: &PairwiseParams) -> DnaAligner {
    let match_fn = MatchParams::new(params.match_score, params.mismatch_score);
    Aligner::with_capacity(
        oligo_len,
        max_ref_len,
        params.gap_open_penalty,
        params.gap_extend_penalty,
        match_fn,
    )
}

#[derive(Debug, Clone)]
pub struct PairwiseMatch {
    pub matched_sequence: String,
    pub score: i32,
    pub mismatches: usize,
    pub has_gaps: bool,
    pub full_coverage: bool,
    /// 0-based start position of the matched window on the reference, when
    /// `has_gaps == false && full_coverage == true`. Zero otherwise.
    pub start: usize,
}

fn process_alignment<F: MatchFunc>(
    aligner: &mut Aligner<F>,
    oligo: &[u8],
    reference: &[u8],
) -> PairwiseMatch {
    let alignment = aligner.local(oligo, reference);

    let mut has_gaps = false;
    let mut mismatches = 0;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {}
            AlignmentOperation::Subst => mismatches += 1,
            AlignmentOperation::Del | AlignmentOperation::Ins => has_gaps = true,
            AlignmentOperation::Xclip(_) | AlignmentOperation::Yclip(_) => {}
        }
    }

    let aligned_query_len = alignment.xend - alignment.xstart;
    let full_coverage = aligned_query_len == oligo.len();

    let matched_sequence = if !has_gaps && full_coverage {
        String::from_utf8_lossy(&reference[alignment.ystart..alignment.yend]).to_string()
    } else {
        String::new()
    };

    PairwiseMatch {
        matched_sequence,
        score: alignment.score,
        mismatches,
        has_gaps,
        full_coverage,
        start: alignment.ystart,
    }
}

pub fn align_oligo_to_reference(
    oligo: &[u8],
    reference: &[u8],
    params: &PairwiseParams,
) -> PairwiseMatch {
    let match_score = params.match_score;
    let mismatch_score = params.mismatch_score;

    let mut aligner = Aligner::with_capacity(
        oligo.len(),
        reference.len(),
        params.gap_open_penalty,
        params.gap_extend_penalty,
        |a: u8, b: u8| -> i32 {
            if a == b { match_score } else { mismatch_score }
        },
    );

    process_alignment(&mut aligner, oligo, reference)
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
        // Anything else (N, IUPAC codes, '-', etc.) is passed through; keeps
        // parity with `simple_align.rs::complement_byte`.
        other => other,
    }
}

fn reverse_complement_bytes(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().rev().map(|&b| complement_byte(b)).collect()
}

fn reverse_complement_str(s: &str) -> String {
    let bytes: Vec<u8> = s.bytes().rev().map(complement_byte).collect();
    String::from_utf8(bytes).unwrap_or_default()
}

/// Run the SW aligner on both the oligo and its reverse complement against
/// the same reference; return the best accepted hit (if any), already
/// re-oriented to the oligo's strand.
///
/// Mirrors the strand-handling semantics of `scan_reference` in
/// `simple_align.rs`:
/// - both strands are tried;
/// - accept criteria are
///   `full_coverage && !has_gaps && mismatches <= max_mismatches`;
/// - tie-break by lowest mismatches, then lowest `ystart` (forward-strand
///   coordinate in the reference), then Forward over Reverse;
/// - when the Reverse pass wins, the matched fragment is reverse-complemented
///   back to oligo orientation.
fn best_strand_match<F: MatchFunc>(
    aligner: &mut Aligner<F>,
    oligo: &[u8],
    rc_oligo: &[u8],
    reference: &[u8],
    max_mismatches: u32,
) -> Option<(PairwiseMatch, Orientation)> {
    let fwd = process_alignment(aligner, oligo, reference);
    let rev = process_alignment(aligner, rc_oligo, reference);

    let fwd_ok = fwd.full_coverage
        && !fwd.has_gaps
        && fwd.mismatches <= max_mismatches as usize;
    let rev_ok = rev.full_coverage
        && !rev.has_gaps
        && rev.mismatches <= max_mismatches as usize;

    match (fwd_ok, rev_ok) {
        (false, false) => None,
        (true, false) => Some((fwd, Orientation::Forward)),
        (false, true) => Some((rc_orient(rev), Orientation::Reverse)),
        (true, true) => {
            let fwd_wins = match fwd.mismatches.cmp(&rev.mismatches) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Greater => false,
                std::cmp::Ordering::Equal => fwd.start <= rev.start,
            };
            if fwd_wins {
                Some((fwd, Orientation::Forward))
            } else {
                Some((rc_orient(rev), Orientation::Reverse))
            }
        }
    }
}

fn rc_orient(mut m: PairwiseMatch) -> PairwiseMatch {
    m.matched_sequence = reverse_complement_str(&m.matched_sequence);
    m
}

pub fn collect_matches(
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> (Vec<String>, usize) {
    let mut matched = Vec::new();
    let mut no_match_count = 0;

    if references.is_empty() {
        return (matched, no_match_count);
    }

    let max_ref_len = references.iter().map(|r| r.len()).max().unwrap();
    let match_score = params.match_score;
    let mismatch_score = params.mismatch_score;

    let mut aligner = Aligner::with_capacity(
        oligo.len(),
        max_ref_len,
        params.gap_open_penalty,
        params.gap_extend_penalty,
        |a: u8, b: u8| -> i32 {
            if a == b { match_score } else { mismatch_score }
        },
    );
    let rc_oligo = reverse_complement_bytes(oligo);

    for reference in references {
        match best_strand_match(&mut aligner, oligo, &rc_oligo, reference, params.max_mismatches) {
            Some((result, _)) => matched.push(result.matched_sequence),
            None => no_match_count += 1,
        }
    }

    (matched, no_match_count)
}

pub fn collect_matches_with_aligner(
    aligner: &mut DnaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> (Vec<String>, usize) {
    let mut matched = Vec::new();
    let mut no_match_count = 0;
    let rc_oligo = reverse_complement_bytes(oligo);

    for reference in references {
        match best_strand_match(aligner, oligo, &rc_oligo, reference, params.max_mismatches) {
            Some((result, _)) => matched.push(result.matched_sequence),
            None => no_match_count += 1,
        }
    }

    (matched, no_match_count)
}

/// Per-reference anchor positions. Same accept rules as
/// `collect_matches_with_aligner`; the returned start is 0-based on the
/// reference's forward strand (matching the bitap path's coordinate system),
/// regardless of which strand the best hit came from. Used by anchored mode.
pub fn collect_anchors_with_aligner(
    aligner: &mut DnaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> Vec<Option<AnchorHit>> {
    let rc_oligo = reverse_complement_bytes(oligo);
    references
        .iter()
        .map(|reference| {
            best_strand_match(aligner, oligo, &rc_oligo, reference, params.max_mismatches).map(
                |(result, orientation)| AnchorHit {
                    start: result.start,
                    orientation,
                    mismatches: result.mismatches as u32,
                },
            )
        })
        .collect()
}

pub fn collect_mismatch_counts_with_aligner(
    aligner: &mut DnaAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &PairwiseParams,
) -> Vec<Option<u32>> {
    let rc_oligo = reverse_complement_bytes(oligo);
    references
        .iter()
        .map(|reference| {
            best_strand_match(aligner, oligo, &rc_oligo, reference, params.max_mismatches)
                .map(|(result, _)| result.mismatches as u32)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> PairwiseParams {
        PairwiseParams::default()
    }

    #[test]
    fn test_exact_match() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_offset() {
        let oligo = b"TATGGTACGT";
        let reference = b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.matched_sequence, "TATGGTACGT");
    }

    #[test]
    fn test_match_with_mismatch() {
        let oligo = b"TATGGTACGT";
        let reference = b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT";
        let result = align_oligo_to_reference(oligo, reference, &default_params());

        assert!(!result.has_gaps);
        assert!(result.full_coverage);
        assert_eq!(result.mismatches, 1);
        assert_eq!(result.matched_sequence, "TATGGTTCGT");
    }

    #[test]
    fn test_collect_matches_from_example() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"AATATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
            b"TATGGTTCGTCATGTTCTAGAAATGGGCTGTTTT".to_vec(),
            b"GTATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec(),
        ];
        let params = default_params();
        let (matched, no_match) = collect_matches(oligo, &references, &params);

        assert_eq!(no_match, 0);
        assert_eq!(matched.len(), 4);
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTACGT").count(), 3);
        assert_eq!(matched.iter().filter(|s| *s == "TATGGTTCGT").count(), 1);
    }

    #[test]
    fn test_max_mismatches_filter() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATGTT".to_vec(),
            b"TATGGTTCGTCATGTT".to_vec(),
        ];
        let mut params = default_params();
        params.max_mismatches = 0;

        let (matched, no_match) = collect_matches(oligo, &references, &params);
        assert_eq!(matched.len(), 1);
        assert_eq!(no_match, 1);
    }

    #[test]
    fn reverse_complement_match_returns_oligo_oriented_fragment() {
        // RC of oligo TATGGTACGT is ACGTACCATA. Embed that on the forward
        // strand of the reference and confirm the returned fragment is RC'd
        // back to the oligo's orientation. Mirrors the bitap test in
        // `simple_align.rs`.
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![b"NNNNACGTACCATANNNN".to_vec()];
        let mut params = default_params();
        params.max_mismatches = 0;

        let (matched, no_match) = collect_matches(oligo, &references, &params);
        assert_eq!(no_match, 0);
        assert_eq!(matched, vec!["TATGGTACGT".to_string()]);
    }

    #[test]
    fn forward_preferred_over_reverse_on_ties() {
        // A palindromic-ish reference where both forward and reverse passes
        // would accept with equal mismatches. Forward must win the tie so
        // the matched fragment comes back unaltered.
        let oligo = b"AAATTT"; // RC is AAATTT — palindrome
        let references: Vec<Vec<u8>> = vec![b"GGAAATTTGG".to_vec()];
        let mut params = default_params();
        params.max_mismatches = 0;

        let (matched, no_match) = collect_matches(oligo, &references, &params);
        assert_eq!(no_match, 0);
        assert_eq!(matched, vec!["AAATTT".to_string()]);
    }

    #[test]
    fn anchor_orientation_recorded_for_reverse_hits() {
        let oligo = b"TATGGTACGT";
        let references: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATG".to_vec(),         // forward exact match
            b"NNNNACGTACCATANNNN".to_vec(),     // reverse exact match
        ];
        let mut params = default_params();
        params.max_mismatches = 0;

        let mut aligner = create_aligner(oligo.len(), 32, &params);
        let anchors = collect_anchors_with_aligner(&mut aligner, oligo, &references, &params);

        assert!(matches!(
            anchors[0],
            Some(AnchorHit { orientation: Orientation::Forward, mismatches: 0, .. })
        ));
        assert!(matches!(
            anchors[1],
            Some(AnchorHit { orientation: Orientation::Reverse, mismatches: 0, .. })
        ));
        // Reverse-pass start is in forward-strand coordinates of the reference.
        assert_eq!(anchors[1].unwrap().start, 4);
    }
}
