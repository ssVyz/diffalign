//! IUPAC alphabet helpers for the search pattern.
//!
//! `iupac_expand` returns the set of {A,C,G,T} bases a single IUPAC code can
//! match, encoded as a 4-bit field. `complement_set` complements such a set
//! (A↔T, C↔G) so we can build a reverse-complement pattern entirely in the
//! 4-bit base-set domain — no string round-trip required.

pub const BIT_A: u8 = 0b0001;
pub const BIT_C: u8 = 0b0010;
pub const BIT_G: u8 = 0b0100;
pub const BIT_T: u8 = 0b1000;

/// Returns Some(set) for a recognized IUPAC byte, None otherwise.
/// `U` collapses to `T` and `I` collapses to `N` per the spec.
pub fn iupac_expand(byte: u8) -> Option<u8> {
    Some(match byte.to_ascii_uppercase() {
        b'A' => BIT_A,
        b'C' => BIT_C,
        b'G' => BIT_G,
        b'T' | b'U' => BIT_T,
        b'R' => BIT_A | BIT_G,
        b'Y' => BIT_C | BIT_T,
        b'S' => BIT_C | BIT_G,
        b'W' => BIT_A | BIT_T,
        b'K' => BIT_G | BIT_T,
        b'M' => BIT_A | BIT_C,
        b'B' => BIT_C | BIT_G | BIT_T,
        b'D' => BIT_A | BIT_G | BIT_T,
        b'H' => BIT_A | BIT_C | BIT_T,
        b'V' => BIT_A | BIT_C | BIT_G,
        b'N' | b'I' => BIT_A | BIT_C | BIT_G | BIT_T,
        _ => return None,
    })
}

/// Complement of a 4-bit base-set: A↔T, C↔G.
pub fn complement_set(set: u8) -> u8 {
    let a = (set & BIT_A) << 3; // A → T
    let c = (set & BIT_C) << 1; // C → G
    let g = (set & BIT_G) >> 1; // G → C
    let t = (set & BIT_T) >> 3; // T → A
    a | c | g | t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_bases() {
        assert_eq!(iupac_expand(b'A'), Some(BIT_A));
        assert_eq!(iupac_expand(b'a'), Some(BIT_A));
        assert_eq!(iupac_expand(b'U'), Some(BIT_T));
        assert_eq!(iupac_expand(b'I'), Some(BIT_A | BIT_C | BIT_G | BIT_T));
        assert_eq!(iupac_expand(b'-'), None);
        assert_eq!(iupac_expand(b'X'), None);
    }

    #[test]
    fn complement_roundtrip() {
        for code in [
            b'A', b'C', b'G', b'T', b'R', b'Y', b'S', b'W', b'K', b'M', b'B',
            b'D', b'H', b'V', b'N',
        ] {
            let set = iupac_expand(code).unwrap();
            assert_eq!(complement_set(complement_set(set)), set, "roundtrip {}", code as char);
        }
        assert_eq!(complement_set(BIT_A), BIT_T);
        assert_eq!(complement_set(BIT_C), BIT_G);
        assert_eq!(complement_set(BIT_A | BIT_G), BIT_T | BIT_C); // R → Y
    }
}
