//! Single-u64 bitap, mismatches-only (Wu–Manber substitution-only variant).
//!
//! ## Boundary contract
//!
//! This module is the simplescreen analogue of `pairwise.rs` in diffalign.
//! The hot path is `scan(...)`. Its inputs are:
//!
//! - `state: &mut BitapState` — per-call reusable scratch (k+1 `u64` registers).
//!   The caller resets it between scans. Holds no other allocations; cheap to
//!   construct, cheap to keep alive across many calls. `Send`, not `Sync`.
//! - `masks: &[u64; 256]` — query-byte → position-bitmask table built once
//!   per pattern by `pattern::PreparedPattern`. For non-ACGT query bytes the
//!   mask is zero, which the recurrence turns into a forced mismatch.
//! - `query: &[u8]` — read-only borrow of one query sequence.
//! - `end_bit: u64` — `1 << (pattern_len - 1)`. The bit we test for a match.
//! - `on_hit: FnMut(end_pos, mismatches)` — called once per text position
//!   where the pattern aligns with ≤ k substitutions, with the smallest
//!   mismatch count for that position.
//!
//! A replacement algorithm in this slot needs to keep the same shape: a
//! per-call worker state, prebuilt pattern data, a `&[u8]` query, and a
//! callback for hits. No allocations on the hot path.
//!
//! ## Algorithm
//!
//! Register `R[d]` is a bitmask where bit `i` is set iff the prefix
//! `p[0..i+1]` matches the text suffix ending at the current position with
//! at most `d` substitutions. Recurrence on text byte `c`:
//!
//!   R[0]_new = ((R[0] << 1) | 1) & B[c]
//!   R[d]_new = ((R[d] << 1) & B[c]) | (R[d-1]_old << 1) | 1   for d ≥ 1
//!
//! Updates run from `d = k` down to `d = 1`, then `d = 0`, so each step still
//! sees the previous step's `R[d-1]` when computing `R[d]`. A match is
//! reported whenever `R[k] & end_bit != 0`; the reported mismatch count is
//! the smallest `d` with `R[d] & end_bit != 0` (since `R[d] ⊆ R[d+1]`).
//!
//! No gaps are allowed by design — the spec is mismatches-only, and the
//! Wu–Manber substitution-only form is strictly faster than the full
//! "k differences" variant.

pub struct BitapState {
    /// `regs[d]` is `R[d]`. Length is `max_mismatches + 1`.
    pub regs: Vec<u64>,
}

impl BitapState {
    pub fn new(max_mismatches: u32) -> Self {
        Self {
            regs: vec![0u64; (max_mismatches + 1) as usize],
        }
    }

    #[inline]
    pub fn reset(&mut self) {
        for r in &mut self.regs {
            *r = 0;
        }
    }
}

/// The hot loop. Walks `query` once, calling `on_hit(end_pos, mismatches)`
/// for every position where the pattern matches with ≤ k substitutions.
///
/// `end_pos` is the 0-based index of the last query byte covered by the
/// match. The match start (0-based) is `end_pos + 1 - pattern_len`.
#[inline]
pub fn scan(
    state: &mut BitapState,
    masks: &[u64; 256],
    query: &[u8],
    end_bit: u64,
    mut on_hit: impl FnMut(u32, u32),
) {
    let regs = state.regs.as_mut_slice();
    let k = regs.len() - 1;

    for (j, &c) in query.iter().enumerate() {
        let cm = masks[c as usize];

        // Update R[k] .. R[1] using R[d-1]_old (still untouched this step).
        for d in (1..=k).rev() {
            regs[d] = ((regs[d] << 1) & cm) | (regs[d - 1] << 1) | 1;
        }
        // R[0] last.
        regs[0] = ((regs[0] << 1) | 1) & cm;

        // Match check: regs[k] is the union (regs[d] ⊆ regs[k] for all d ≤ k).
        if regs[k] & end_bit != 0 {
            // Find the smallest d with the end bit set — that's the actual
            // mismatch count of this hit.
            for d in 0..=k {
                if regs[d] & end_bit != 0 {
                    on_hit(j as u32, d as u32);
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::pattern::PreparedPattern;

    fn collect_fwd(pattern: &[u8], k: u32, query: &[u8]) -> Vec<(u32, u32)> {
        let p = PreparedPattern::build(pattern, k).unwrap();
        let mut s = BitapState::new(k);
        let mut hits = Vec::new();
        scan(&mut s, &p.fwd_masks, query, p.end_bit, |end, mm| {
            hits.push((end, mm));
        });
        hits
    }

    #[test]
    fn exact_match() {
        let hits = collect_fwd(b"ACGT", 0, b"AAACGTAA");
        assert_eq!(hits, vec![(5, 0)]); // end_pos=5 → start=2
    }

    #[test]
    fn one_mismatch_finds_zero_and_one() {
        let hits = collect_fwd(b"ACGT", 1, b"ACGTACAT");
        // end=3: ACGT vs ACGT — 0 mm
        // end=7: ACAT vs ACGT — 1 mm
        // end=4..6 may not match (frame-shifted alignments produce many
        // mismatches very fast).
        assert!(hits.contains(&(3, 0)));
        assert!(hits.contains(&(7, 1)));
    }

    #[test]
    fn no_hit_when_too_many_mismatches() {
        let hits = collect_fwd(b"ACGT", 1, b"GGGG");
        assert!(hits.is_empty());
    }

    #[test]
    fn shorter_query_than_pattern() {
        let hits = collect_fwd(b"ACGTACGT", 0, b"ACG");
        assert!(hits.is_empty());
    }

    #[test]
    fn iupac_pattern_matches_multiple_bases() {
        // Pattern "RCGT" (R=A|G). Should match both ACGT and GCGT exactly.
        let hits = collect_fwd(b"RCGT", 0, b"ACGTNNNGCGT");
        assert!(hits.contains(&(3, 0)));
        assert!(hits.contains(&(10, 0)));
    }

    #[test]
    fn n_in_query_is_a_mismatch() {
        // 0 mismatches allowed: a single N inside the window blocks the match.
        let hits = collect_fwd(b"ACGT", 0, b"ACNT");
        assert!(hits.is_empty());
        // 1 mismatch allowed: the N counts as one.
        let hits = collect_fwd(b"ACGT", 1, b"ACNT");
        assert!(hits.contains(&(3, 1)));
    }

    #[test]
    fn reports_smallest_mismatch_count() {
        // Pattern matches at end=3 with 0 mm; with k=2 the same end position
        // is still reported as mm=0, not mm=2.
        let hits = collect_fwd(b"ACGT", 2, b"ACGT");
        assert_eq!(hits, vec![(3, 0)]);
    }

    #[test]
    fn invalid_pattern_yields_no_hits() {
        let p = PreparedPattern::build(b"ACX T", 0).unwrap();
        assert!(!p.valid);
        let mut s = BitapState::new(0);
        let mut hits = Vec::new();
        scan(&mut s, &p.fwd_masks, b"ACGTACGT", p.end_bit, |e, m| hits.push((e, m)));
        assert!(hits.is_empty());
    }

    #[test]
    fn pattern_at_max_length_works() {
        let pat = vec![b'A'; 64];
        let mut q = vec![b'C'; 100];
        for i in 10..74 {
            q[i] = b'A';
        }
        let hits = collect_fwd(&pat, 0, &q);
        assert_eq!(hits, vec![(73, 0)]); // start = 73 - 64 + 1 = 10
    }
}
