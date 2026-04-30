//! Per-query orchestration: drive the bitap forward, drive it on the
//! reverse-complement masks, return all hits with 1-based forward-strand
//! coordinates.

use super::bitap::{BitapState, scan};
use super::pattern::PreparedPattern;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    pub fn symbol(self) -> char {
        match self {
            Orientation::Forward => '+',
            Orientation::Reverse => '-',
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Hit {
    /// 1-based start coordinate on the forward strand of the query (inclusive).
    pub start: u32,
    /// 1-based end coordinate on the forward strand of the query (inclusive).
    pub end: u32,
    pub mismatches: u32,
    pub orientation: Orientation,
}

/// Run both directions over `query` and append hits to `out`. The caller
/// owns `out` and may pre-allocate / reuse it across many queries.
///
/// `state` is per-call scratch — same `BitapState` may be passed for every
/// query in a run. It is reset internally between the forward and reverse
/// passes.
pub fn screen(
    state: &mut BitapState,
    pattern: &PreparedPattern,
    query: &[u8],
    out: &mut Vec<Hit>,
) {
    if !pattern.valid {
        return;
    }

    let len = pattern.len;
    let end_bit = pattern.end_bit;

    state.reset();
    scan(state, &pattern.fwd_masks, query, end_bit, |end_pos, mm| {
        // end_pos is 0-based. The end-bit cannot fire before `len` bytes have
        // been consumed, so end_pos + 1 >= len always.
        let start_0 = end_pos + 1 - len;
        out.push(Hit {
            start: start_0 + 1,
            end: end_pos + 1,
            mismatches: mm,
            orientation: Orientation::Forward,
        });
    });

    state.reset();
    scan(state, &pattern.rc_masks, query, end_bit, |end_pos, mm| {
        let start_0 = end_pos + 1 - len;
        out.push(Hit {
            start: start_0 + 1,
            end: end_pos + 1,
            mismatches: mm,
            orientation: Orientation::Reverse,
        });
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_and_reverse_complement_both_found() {
        // Pattern: ACGT. RC: ACGT (palindrome). Query has it once.
        let p = PreparedPattern::build(b"ACGT", 0).unwrap();
        let mut s = BitapState::new(0);
        let mut hits = Vec::new();
        screen(&mut s, &p, b"NNNACGTNNN", &mut hits);
        // Both fwd and rev should fire at the same site (palindrome).
        assert_eq!(hits.len(), 2);
        assert!(hits.iter().any(|h| h.orientation == Orientation::Forward && h.start == 4));
        assert!(hits.iter().any(|h| h.orientation == Orientation::Reverse && h.start == 4));
    }

    #[test]
    fn reverse_complement_found_on_correct_strand_coordinate() {
        // Pattern: ACGG. RC: CCGT.
        // Place CCGT at query positions 4..7 (1-based 5..8).
        let p = PreparedPattern::build(b"ACGG", 0).unwrap();
        let mut s = BitapState::new(0);
        let mut hits = Vec::new();
        screen(&mut s, &p, b"NNNNCCGTNN", &mut hits);
        assert_eq!(hits.len(), 1);
        let h = hits[0];
        assert_eq!(h.orientation, Orientation::Reverse);
        assert_eq!(h.start, 5);
        assert_eq!(h.end, 8);
        assert_eq!(h.mismatches, 0);
    }

    #[test]
    fn invalid_pattern_yields_no_hits() {
        let p = PreparedPattern::build(b"ACX T", 0).unwrap();
        let mut s = BitapState::new(0);
        let mut hits = Vec::new();
        screen(&mut s, &p, b"ACGTACGT", &mut hits);
        assert!(hits.is_empty());
    }
}
