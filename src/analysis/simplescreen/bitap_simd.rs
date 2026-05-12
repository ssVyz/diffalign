//! AVX2-vectorized bitap. Processes `LANES_AVX2` references concurrently with
//! the same pattern.
//!
//! Algorithm and correctness invariants are identical to the scalar `bitap.rs`:
//! the recurrence on register `R[d]` is performed in parallel across 4 lanes
//! (each lane holds one reference's state in a 64-bit slot of an `__m256i`).
//!
//! ## Why this matches scalar output bit-for-bit
//!
//! - The bitap recurrence is exact integer arithmetic on bitmasks — no
//!   nondeterminism, no rounding.
//! - Each lane's state evolves only from its own previous state and its own
//!   per-step byte; lanes never interact.
//! - Hits are tracked per lane with the same `(min mismatches, earliest end
//!   position)` rule the scalar code's downstream `best_hit` selection would
//!   produce for a single direction.
//!
//! ## Ragged-length handling
//!
//! References within a 4-batch typically have different lengths. We iterate
//! `j` from 0 to `max(lens[*])`; once a lane has run past its reference end we
//! feed byte `0` to that lane (mask becomes 0, forcing mismatches). That state
//! can in principle re-trigger `end_bit` from leftover register bits, so we
//! also suppress hit *recording* for any lane whose `j >= lens[lane]`. With
//! both guards in place, padded steps never contribute spurious hits.

#![cfg(target_arch = "x86_64")]

use std::arch::x86_64::*;

pub const LANES_AVX2: usize = 4;

/// Pattern length cap is 64 bp, so `max_mismatches < 64` and we need at most
/// 64 registers. Allocated on the stack at this fixed size.
const MAX_REGS: usize = 64;

/// Per-lane result of a single-direction scan.
#[derive(Copy, Clone, Debug, Default)]
pub struct LaneBest {
    pub mismatches: u32,
    /// 0-based end position of the best hit found in this lane.
    pub end_pos: u32,
    pub found: bool,
}

/// Scan a 4-tuple of references with one mask table, returning the best hit
/// per lane (minimum mismatches; ties broken by earliest end position).
///
/// # Safety
/// Caller must have verified at runtime that the CPU supports AVX2.
#[target_feature(enable = "avx2")]
pub unsafe fn scan_avx2(
    masks: &[u64; 256],
    refs: &[&[u8]; LANES_AVX2],
    max_mismatches: u32,
    end_bit: u64,
) -> [LaneBest; LANES_AVX2] {
    let k = max_mismatches as usize;
    debug_assert!(k + 1 <= MAX_REGS);

    let mut regs = [_mm256_setzero_si256(); MAX_REGS];
    let one_v = _mm256_set1_epi64x(1);
    let end_bit_v = _mm256_set1_epi64x(end_bit as i64);

    let lens: [usize; LANES_AVX2] = [refs[0].len(), refs[1].len(), refs[2].len(), refs[3].len()];
    let max_len = lens.iter().copied().max().unwrap_or(0);

    let mut best = [LaneBest::default(); LANES_AVX2];

    let ptrs: [*const u8; LANES_AVX2] = [
        refs[0].as_ptr(),
        refs[1].as_ptr(),
        refs[2].as_ptr(),
        refs[3].as_ptr(),
    ];

    for j in 0..max_len {
        let b0 = if j < lens[0] { *ptrs[0].add(j) } else { 0 };
        let b1 = if j < lens[1] { *ptrs[1].add(j) } else { 0 };
        let b2 = if j < lens[2] { *ptrs[2].add(j) } else { 0 };
        let b3 = if j < lens[3] { *ptrs[3].add(j) } else { 0 };

        // `_mm256_set_epi64x` packs in reverse: high lane first.
        let cm = _mm256_set_epi64x(
            masks[b3 as usize] as i64,
            masks[b2 as usize] as i64,
            masks[b1 as usize] as i64,
            masks[b0 as usize] as i64,
        );

        // R[d]_new = ((R[d] << 1) & cm) | (R[d-1]_old << 1) | 1, for d >= 1.
        // Iterate from d = k down to d = 1 so each step still sees the
        // previous step's `R[d-1]` (mirrors scalar bitap ordering).
        for d in (1..=k).rev() {
            let shifted_d = _mm256_slli_epi64::<1>(regs[d]);
            let shifted_dm1 = _mm256_slli_epi64::<1>(regs[d - 1]);
            let a = _mm256_and_si256(shifted_d, cm);
            let b = _mm256_or_si256(a, shifted_dm1);
            regs[d] = _mm256_or_si256(b, one_v);
        }
        // R[0]_new = ((R[0] << 1) | 1) & cm.
        let shifted_0 = _mm256_slli_epi64::<1>(regs[0]);
        let with_one = _mm256_or_si256(shifted_0, one_v);
        regs[0] = _mm256_and_si256(with_one, cm);

        // Cheap "any lane hit" check via VPTEST.
        let hits_k = _mm256_and_si256(regs[k], end_bit_v);
        if _mm256_testz_si256(hits_k, hits_k) != 0 {
            continue;
        }

        // At least one lane hit at d=k. Extract regs[0..=k] to find the
        // smallest d per hitting lane.
        let mut reg_arrs = [[0u64; LANES_AVX2]; MAX_REGS];
        for d in 0..=k {
            _mm256_storeu_si256(reg_arrs[d].as_mut_ptr() as *mut __m256i, regs[d]);
        }

        for lane in 0..LANES_AVX2 {
            if j >= lens[lane] {
                continue;
            }
            if reg_arrs[k][lane] & end_bit == 0 {
                continue;
            }
            // Once a lane has a 0-mm hit at its earliest position, it can't
            // improve further; skip the per-d search.
            if best[lane].found && best[lane].mismatches == 0 {
                continue;
            }
            for d in 0..=k {
                if reg_arrs[d][lane] & end_bit != 0 {
                    let d_u32 = d as u32;
                    if !best[lane].found || d_u32 < best[lane].mismatches {
                        best[lane] = LaneBest {
                            mismatches: d_u32,
                            end_pos: j as u32,
                            found: true,
                        };
                    }
                    // Smallest d wins; no need to check higher d.
                    break;
                }
            }
        }
    }

    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::simplescreen::pattern::PreparedPattern;

    fn skip_if_no_avx2() -> bool {
        if std::is_x86_feature_detected!("avx2") {
            false
        } else {
            eprintln!("skipping AVX2 bitap test on a CPU without AVX2");
            true
        }
    }

    #[test]
    fn scan_avx2_matches_scalar_one_lane() {
        if skip_if_no_avx2() {
            return;
        }
        // Lane 0 has a single exact match at end_pos = 9. Other lanes empty.
        let pattern = PreparedPattern::build(b"TATGGTACGT", 1).unwrap();
        let lane0: &[u8] = b"TATGGTACGTAAAA";
        let empty: &[u8] = b"";
        let refs: [&[u8]; LANES_AVX2] = [lane0, empty, empty, empty];
        let out = unsafe { scan_avx2(&pattern.fwd_masks, &refs, 1, pattern.end_bit) };
        assert!(out[0].found);
        assert_eq!(out[0].mismatches, 0);
        assert_eq!(out[0].end_pos, 9);
        assert!(!out[1].found);
        assert!(!out[2].found);
        assert!(!out[3].found);
    }

    #[test]
    fn scan_avx2_four_lanes_independent_results() {
        if skip_if_no_avx2() {
            return;
        }
        let pattern = PreparedPattern::build(b"ACGT", 1).unwrap();
        let l0: &[u8] = b"NNACGTNN";        // exact at end=5
        let l1: &[u8] = b"NNACCTNN";        // 1mm at end=5
        let l2: &[u8] = b"NNNNNNNN";        // no match
        let l3: &[u8] = b"NNACGTACGT";      // two exact; earliest is end=5
        let refs: [&[u8]; LANES_AVX2] = [l0, l1, l2, l3];
        let out = unsafe { scan_avx2(&pattern.fwd_masks, &refs, 1, pattern.end_bit) };
        assert!(out[0].found);
        assert_eq!((out[0].mismatches, out[0].end_pos), (0, 5));
        assert!(out[1].found);
        assert_eq!((out[1].mismatches, out[1].end_pos), (1, 5));
        assert!(!out[2].found);
        assert!(out[3].found);
        assert_eq!((out[3].mismatches, out[3].end_pos), (0, 5));
    }

    #[test]
    fn scan_avx2_picks_lowest_mismatch_then_earliest() {
        if skip_if_no_avx2() {
            return;
        }
        let pattern = PreparedPattern::build(b"ACGT", 2).unwrap();
        // Two hits in this lane: 1-mm at end=3, exact at end=10.
        // Best by tiebreaker is the exact-match (lowest mm wins).
        let r0: &[u8] = b"ACATXXXACGT";
        let empty: &[u8] = b"";
        let refs: [&[u8]; LANES_AVX2] = [r0, empty, empty, empty];
        let out = unsafe { scan_avx2(&pattern.fwd_masks, &refs, 2, pattern.end_bit) };
        assert!(out[0].found);
        assert_eq!(out[0].mismatches, 0);
        assert_eq!(out[0].end_pos, 10);
    }

    #[test]
    fn scan_avx2_padding_does_not_emit_spurious_hits() {
        if skip_if_no_avx2() {
            return;
        }
        // Reference is shorter than pattern in every lane → no hit possible.
        let pattern = PreparedPattern::build(b"ACGTACGT", 2).unwrap();
        let short: &[u8] = b"ACGT";
        let refs: [&[u8]; LANES_AVX2] = [short, short, short, short];
        let out = unsafe { scan_avx2(&pattern.fwd_masks, &refs, 2, pattern.end_bit) };
        for r in out.iter() {
            assert!(!r.found);
        }
    }
}
