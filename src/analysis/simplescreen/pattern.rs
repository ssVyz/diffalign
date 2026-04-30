//! Bitap masks pre-built from the search sequence (forward + reverse complement).
//!
//! The hot loop sees only `PreparedPattern` and `bitap::BitapState`. All
//! IUPAC interpretation, validation, and reverse-complement work happens here,
//! exactly once per run.

use super::iupac::{BIT_A, BIT_C, BIT_G, BIT_T, complement_set, iupac_expand};

pub const MAX_PATTERN_LEN: usize = 64;

pub struct PreparedPattern {
    pub len: u32,
    /// Kept on the struct so a downstream caller (e.g. diffalign) can read
    /// back the `k` the masks were built for without tracking it separately.
    #[allow(dead_code)]
    pub max_mismatches: u32,
    /// False if the search sequence contained any byte we couldn't interpret.
    /// In that case both mask tables are all zero and no hits will fire.
    pub valid: bool,
    /// `fwd_masks[c]` has bit i set iff the pattern's i-th position can match
    /// query byte `c`. Lowercase ACGT slots are populated alongside uppercase
    /// so the hot loop never needs to case-fold.
    pub fwd_masks: [u64; 256],
    pub rc_masks: [u64; 256],
    /// `1 << (len - 1)` — the bit we test in each register for "match found".
    pub end_bit: u64,
}

impl PreparedPattern {
    pub fn build(pattern: &[u8], max_mismatches: u32) -> Result<Self, String> {
        let len = pattern.len();
        if len == 0 {
            return Err("search sequence is empty".into());
        }
        if len > MAX_PATTERN_LEN {
            return Err(format!(
                "search sequence is {} bp; this build supports up to {} bp (single-u64 bitap)",
                len, MAX_PATTERN_LEN
            ));
        }
        if max_mismatches as usize >= len {
            return Err(format!(
                "max-mismatches ({}) must be less than the search sequence length ({})",
                max_mismatches, len
            ));
        }

        let mut fwd_set = [0u8; MAX_PATTERN_LEN];
        let mut valid = true;
        for (i, &b) in pattern.iter().enumerate() {
            match iupac_expand(b) {
                Some(s) => fwd_set[i] = s,
                None => {
                    valid = false;
                    break;
                }
            }
        }

        let mut fwd_masks = [0u64; 256];
        let mut rc_masks = [0u64; 256];

        if valid {
            // RC pattern: reverse and complement the per-position base sets.
            let mut rc_set = [0u8; MAX_PATTERN_LEN];
            for i in 0..len {
                rc_set[i] = complement_set(fwd_set[len - 1 - i]);
            }
            fill_masks(&fwd_set[..len], &mut fwd_masks);
            fill_masks(&rc_set[..len], &mut rc_masks);
        }

        Ok(Self {
            len: len as u32,
            max_mismatches,
            valid,
            fwd_masks,
            rc_masks,
            end_bit: 1u64 << (len - 1),
        })
    }
}

fn fill_masks(position_sets: &[u8], masks: &mut [u64; 256]) {
    for (i, &set) in position_sets.iter().enumerate() {
        let bit = 1u64 << i;
        if set & BIT_A != 0 {
            masks[b'A' as usize] |= bit;
            masks[b'a' as usize] |= bit;
        }
        if set & BIT_C != 0 {
            masks[b'C' as usize] |= bit;
            masks[b'c' as usize] |= bit;
        }
        if set & BIT_G != 0 {
            masks[b'G' as usize] |= bit;
            masks[b'g' as usize] |= bit;
        }
        if set & BIT_T != 0 {
            masks[b'T' as usize] |= bit;
            masks[b't' as usize] |= bit;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwd_masks_are_position_indicators() {
        let p = PreparedPattern::build(b"ACGT", 0).unwrap();
        assert!(p.valid);
        assert_eq!(p.len, 4);
        assert_eq!(p.end_bit, 0b1000);
        assert_eq!(p.fwd_masks[b'A' as usize], 0b0001);
        assert_eq!(p.fwd_masks[b'C' as usize], 0b0010);
        assert_eq!(p.fwd_masks[b'G' as usize], 0b0100);
        assert_eq!(p.fwd_masks[b'T' as usize], 0b1000);
        assert_eq!(p.fwd_masks[b'N' as usize], 0); // N in query: zero mask.
    }

    #[test]
    fn rc_masks_match_reverse_complement() {
        // pattern ACGT, RC = ACGT (palindrome), so RC masks equal fwd masks.
        let p = PreparedPattern::build(b"ACGT", 0).unwrap();
        for c in 0u8..=255 {
            assert_eq!(p.fwd_masks[c as usize], p.rc_masks[c as usize], "byte {}", c);
        }
        // pattern AAAA, RC = TTTT.
        let p = PreparedPattern::build(b"AAAA", 0).unwrap();
        assert_eq!(p.fwd_masks[b'A' as usize], 0b1111);
        assert_eq!(p.fwd_masks[b'T' as usize], 0);
        assert_eq!(p.rc_masks[b'A' as usize], 0);
        assert_eq!(p.rc_masks[b'T' as usize], 0b1111);
    }

    #[test]
    fn iupac_expansion_in_pattern() {
        // pattern "RY" (R=A|G, Y=C|T) at positions 0 and 1.
        let p = PreparedPattern::build(b"RY", 0).unwrap();
        assert_eq!(p.fwd_masks[b'A' as usize], 0b01); // matches position 0 only
        assert_eq!(p.fwd_masks[b'G' as usize], 0b01);
        assert_eq!(p.fwd_masks[b'C' as usize], 0b10); // matches position 1 only
        assert_eq!(p.fwd_masks[b'T' as usize], 0b10);
    }

    #[test]
    fn lowercase_query_bytes_are_handled() {
        let p = PreparedPattern::build(b"ACGT", 0).unwrap();
        assert_eq!(p.fwd_masks[b'a' as usize], p.fwd_masks[b'A' as usize]);
        assert_eq!(p.fwd_masks[b'g' as usize], p.fwd_masks[b'G' as usize]);
    }

    #[test]
    fn unknown_byte_in_pattern_invalidates() {
        let p = PreparedPattern::build(b"ACX T", 0).unwrap();
        assert!(!p.valid);
        // All masks zero — no hits will fire.
        for c in 0u8..=255 {
            assert_eq!(p.fwd_masks[c as usize], 0);
            assert_eq!(p.rc_masks[c as usize], 0);
        }
    }

    #[test]
    fn rejects_too_long() {
        let pat = vec![b'A'; 65];
        assert!(PreparedPattern::build(&pat, 0).is_err());
    }

    #[test]
    fn rejects_max_mismatches_at_or_above_len() {
        assert!(PreparedPattern::build(b"ACGT", 4).is_err());
        assert!(PreparedPattern::build(b"ACGT", 5).is_err());
        assert!(PreparedPattern::build(b"ACGT", 3).is_ok());
    }
}
