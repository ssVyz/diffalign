//! Bitap-based aligner backend (simplescreen).
//!
//! Mirrors the API surface of `pairwise.rs`:
//!
//! * `SimpleAligner` — per-worker reusable state (analogue of `DnaAligner`).
//! * `create_simple_aligner` — constructed once via Rayon's `map_init`.
//! * `collect_matches_with_simple_aligner` — oligo vs all references, returns
//!   matched fragments + no-match count.
//! * `collect_mismatch_counts_with_simple_aligner` — per-reference mismatch
//!   count (`None` = rejected).
//!
//! Forward + reverse-complement scanning is enabled. When the best hit for a
//! reference is on the reverse strand, the matched fragment is
//! reverse-complemented before being returned, so downstream variant analysis
//! sees all fragments in the same orientation as the oligo.

use super::simplescreen::bitap::BitapState;
use super::simplescreen::pattern::{MAX_PATTERN_LEN, PreparedPattern};
use super::simplescreen::screener::{Hit, Orientation, screen};
use super::types::{AnchorHit, SimpleParams};

/// Per-worker scratch state for the bitap aligner.
///
/// `Send` (the inner `Vec<u64>` and `Vec<Hit>` are `Send`); not `Sync`, which
/// is fine — Rayon hands one `SimpleAligner` to each worker via `map_init`.
pub struct SimpleAligner {
    state: BitapState,
    /// Reusable hit buffer; cleared per reference scan.
    hits: Vec<Hit>,
}

pub fn create_simple_aligner(params: &SimpleParams) -> SimpleAligner {
    SimpleAligner {
        state: BitapState::new(params.max_mismatches),
        hits: Vec::with_capacity(16),
    }
}

/// Largest oligo length the bitap backend can handle.
pub const SIMPLE_MAX_OLIGO_LEN: usize = MAX_PATTERN_LEN;

#[derive(Debug, Clone)]
pub struct SimpleMatch {
    /// The matched fragment, oriented to match the oligo's strand
    /// (reverse-complemented if the best hit was on the reverse strand).
    pub matched_sequence: String,
    pub mismatches: u32,
    pub orientation: Orientation,
}

/// Pick the best hit from the buffer: lowest mismatches, then lowest start.
/// Forward orientation breaks remaining ties (deterministic, mirrors SW's
/// "leftmost-best" tendency).
fn best_hit(hits: &[Hit]) -> Option<Hit> {
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
        // Anything else (N, IUPAC codes, '-', etc.) is passed through; these
        // only show up in references and won't affect grouping correctness.
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

/// Run forward+RC scan on one reference and return the best hit (if any),
/// already extracting the matched substring in oligo orientation.
fn scan_reference(
    aligner: &mut SimpleAligner,
    pattern: &PreparedPattern,
    reference: &[u8],
) -> Option<SimpleMatch> {
    aligner.hits.clear();
    screen(&mut aligner.state, pattern, reference, &mut aligner.hits);

    let hit = best_hit(&aligner.hits)?;
    let len = pattern.len as usize;
    // `start` and `end` are 1-based, inclusive; convert to a [start_0..end_0)
    // half-open slice on the forward strand.
    let start_0 = (hit.start as usize) - 1;
    let end_0 = hit.end as usize;
    let frag = &reference[start_0..end_0];
    debug_assert_eq!(frag.len(), len);

    let matched_sequence = match hit.orientation {
        Orientation::Forward => String::from_utf8_lossy(frag).into_owned(),
        Orientation::Reverse => reverse_complement(frag),
    };

    Some(SimpleMatch {
        matched_sequence,
        mismatches: hit.mismatches,
        orientation: hit.orientation,
    })
}

pub fn collect_matches_with_simple_aligner(
    aligner: &mut SimpleAligner,
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

    let mut matched = Vec::new();
    let mut no_match_count = 0usize;
    for reference in references {
        match scan_reference(aligner, &pattern, reference) {
            Some(m) => matched.push(m.matched_sequence),
            None => no_match_count += 1,
        }
    }
    (matched, no_match_count)
}

/// Per-reference anchor positions. Same accept rules as
/// `collect_matches_with_simple_aligner`; the returned start position is
/// 0-based on the reference's forward strand (same coordinate system the
/// downstream length-`L` fragment extractor expects).
pub fn collect_anchors_with_simple_aligner(
    aligner: &mut SimpleAligner,
    oligo: &[u8],
    references: &[Vec<u8>],
    params: &SimpleParams,
) -> Vec<Option<AnchorHit>> {
    let pattern = match PreparedPattern::build(oligo, params.max_mismatches) {
        Ok(p) => p,
        Err(_) => return references.iter().map(|_| None).collect(),
    };
    if !pattern.valid {
        return references.iter().map(|_| None).collect();
    }

    references
        .iter()
        .map(|reference| {
            aligner.hits.clear();
            screen(&mut aligner.state, &pattern, reference, &mut aligner.hits);
            best_hit(&aligner.hits).map(|h| AnchorHit {
                start: (h.start as usize) - 1,
                orientation: h.orientation,
                mismatches: h.mismatches,
            })
        })
        .collect()
}

pub fn collect_mismatch_counts_with_simple_aligner(
    aligner: &mut SimpleAligner,
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

    references
        .iter()
        .map(|reference| scan_reference(aligner, &pattern, reference).map(|m| m.mismatches))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> SimpleParams {
        SimpleParams { max_mismatches: 2 }
    }

    #[test]
    fn forward_exact_match() {
        let mut aligner = create_simple_aligner(&p());
        let oligo = b"TATGGTACGT";
        let refs: Vec<Vec<u8>> = vec![b"TATGGTACGTCATGTTCTAGAAATGGGCTGT".to_vec()];
        let (matched, nm) = collect_matches_with_simple_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 0);
        assert_eq!(matched, vec!["TATGGTACGT".to_string()]);
    }

    #[test]
    fn reverse_complement_match_returns_oligo_oriented_fragment() {
        let mut aligner = create_simple_aligner(&p());
        let oligo = b"TATGGTACGT";
        // RC of oligo = ACGTACCATA. Place that on the forward strand.
        let refs: Vec<Vec<u8>> = vec![b"NNNNACGTACCATANNNN".to_vec()];
        let (matched, nm) = collect_matches_with_simple_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 0);
        // The returned fragment should be RC'd back to the oligo's orientation.
        assert_eq!(matched, vec!["TATGGTACGT".to_string()]);
    }

    #[test]
    fn picks_best_when_multiple_hits() {
        let mut aligner = create_simple_aligner(&p());
        let oligo = b"ACGT";
        // Two forward placements: at pos 0 (mm=1, ACGA→ACGT differs by 1)
        // and at pos 5 (mm=0, ACGT). Pick the 0-mm one.
        let refs: Vec<Vec<u8>> = vec![b"ACGAANACGTNN".to_vec()];
        let mut params = p();
        params.max_mismatches = 1;
        let (matched, nm) = collect_matches_with_simple_aligner(&mut aligner, oligo, &refs, &params);
        assert_eq!(nm, 0);
        assert_eq!(matched, vec!["ACGT".to_string()]);
    }

    #[test]
    fn rejects_above_max_mismatches() {
        let mut aligner = create_simple_aligner(&p());
        let oligo = b"ACGTACGTACGT";
        let refs: Vec<Vec<u8>> = vec![b"GGGGGGGGGGGG".to_vec()];
        let (matched, nm) = collect_matches_with_simple_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(nm, 1);
        assert!(matched.is_empty());
    }

    #[test]
    fn mismatch_counts_sized_to_input() {
        let mut aligner = create_simple_aligner(&p());
        let oligo = b"TATGGTACGT";
        let refs: Vec<Vec<u8>> = vec![
            b"TATGGTACGTCATG".to_vec(), // exact
            b"GGGGGGGGGGGGGG".to_vec(), // no match
        ];
        let counts =
            collect_mismatch_counts_with_simple_aligner(&mut aligner, oligo, &refs, &p());
        assert_eq!(counts.len(), 2);
        assert_eq!(counts[0], Some(0));
        assert_eq!(counts[1], None);
    }
}
