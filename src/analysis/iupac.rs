//! IUPAC ambiguity codes and DNA sequence utilities

use once_cell::sync::Lazy;
use std::collections::{HashMap, HashSet};

pub const STANDARD_BASES: [char; 4] = ['A', 'C', 'G', 'T'];

pub const AMBIGUOUS_BASES: [char; 11] = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'];

pub const GAP_CHARS: [char; 2] = ['-', '.'];

pub static IUPAC_TO_BASES: Lazy<HashMap<char, HashSet<char>>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', ['A'].into_iter().collect());
    map.insert('C', ['C'].into_iter().collect());
    map.insert('G', ['G'].into_iter().collect());
    map.insert('T', ['T'].into_iter().collect());
    map.insert('R', ['A', 'G'].into_iter().collect());
    map.insert('Y', ['C', 'T'].into_iter().collect());
    map.insert('S', ['G', 'C'].into_iter().collect());
    map.insert('W', ['A', 'T'].into_iter().collect());
    map.insert('K', ['G', 'T'].into_iter().collect());
    map.insert('M', ['A', 'C'].into_iter().collect());
    map.insert('B', ['C', 'G', 'T'].into_iter().collect());
    map.insert('D', ['A', 'G', 'T'].into_iter().collect());
    map.insert('H', ['A', 'C', 'T'].into_iter().collect());
    map.insert('V', ['A', 'C', 'G'].into_iter().collect());
    map.insert('N', ['A', 'C', 'G', 'T'].into_iter().collect());
    map
});

pub static BASES_TO_IUPAC: Lazy<HashMap<Vec<char>, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert(vec!['A'], 'A');
    map.insert(vec!['C'], 'C');
    map.insert(vec!['G'], 'G');
    map.insert(vec!['T'], 'T');
    map.insert(vec!['A', 'G'], 'R');
    map.insert(vec!['C', 'T'], 'Y');
    map.insert(vec!['C', 'G'], 'S');
    map.insert(vec!['A', 'T'], 'W');
    map.insert(vec!['G', 'T'], 'K');
    map.insert(vec!['A', 'C'], 'M');
    map.insert(vec!['C', 'G', 'T'], 'B');
    map.insert(vec!['A', 'G', 'T'], 'D');
    map.insert(vec!['A', 'C', 'T'], 'H');
    map.insert(vec!['A', 'C', 'G'], 'V');
    map.insert(vec!['A', 'C', 'G', 'T'], 'N');
    map
});

pub static COMPLEMENT: Lazy<HashMap<char, char>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert('A', 'T');
    map.insert('T', 'A');
    map.insert('C', 'G');
    map.insert('G', 'C');
    map.insert('R', 'Y');
    map.insert('Y', 'R');
    map.insert('S', 'S');
    map.insert('W', 'W');
    map.insert('K', 'M');
    map.insert('M', 'K');
    map.insert('B', 'V');
    map.insert('V', 'B');
    map.insert('D', 'H');
    map.insert('H', 'D');
    map.insert('N', 'N');
    map
});

pub fn is_standard_base(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

pub fn is_ambiguous_base(c: char) -> bool {
    matches!(c, 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N')
}

pub fn is_gap(c: char) -> bool {
    matches!(c, '-' | '.')
}

pub fn is_valid_dna(c: char) -> bool {
    is_standard_base(c) || is_ambiguous_base(c)
}

pub fn bases_to_iupac(bases: &HashSet<char>) -> char {
    let mut sorted: Vec<char> = bases.iter().copied().collect();
    sorted.sort();
    *BASES_TO_IUPAC.get(&sorted).unwrap_or(&'N')
}

pub fn iupac_to_bases(code: char) -> Option<&'static HashSet<char>> {
    IUPAC_TO_BASES.get(&code)
}

pub fn sequence_matches_consensus(seq: &str, consensus: &str) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }

    for (s, c) in seq.chars().zip(consensus.chars()) {
        if let Some(allowed) = IUPAC_TO_BASES.get(&c) {
            if !allowed.contains(&s) {
                return false;
            }
        } else if s != c {
            return false;
        }
    }
    true
}

pub fn create_consensus(sequences: &[&str], exclude_n: bool) -> (String, usize, bool) {
    if sequences.is_empty() {
        return (String::new(), 0, true);
    }

    let seq_len = sequences[0].len();
    let mut consensus = String::with_capacity(seq_len);
    let mut ambiguity_count = 0;

    for pos in 0..seq_len {
        let mut bases_at_pos: HashSet<char> = HashSet::new();
        for seq in sequences {
            if let Some(c) = seq.chars().nth(pos) {
                bases_at_pos.insert(c);
            }
        }

        if bases_at_pos.len() == 1 {
            consensus.push(*bases_at_pos.iter().next().unwrap());
        } else {
            let code = bases_to_iupac(&bases_at_pos);
            if exclude_n && code == 'N' {
                return (consensus, ambiguity_count, false);
            }
            consensus.push(code);
            ambiguity_count += 1;
        }
    }

    (consensus, ambiguity_count, true)
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| *COMPLEMENT.get(&c).unwrap_or(&c))
        .collect()
}

pub fn count_ambiguities(seq: &str) -> usize {
    seq.chars().filter(|&c| is_ambiguous_base(c)).count()
}

pub const IUPAC_FROM_MASK: [u8; 16] = [
    b'?', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

#[inline]
pub fn base_to_bit(b: u8) -> u8 {
    match b {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

#[inline]
pub fn iupac_to_mask(b: u8) -> u8 {
    base_to_bit(b)
}

#[inline]
pub fn sequence_matches_consensus_bytes(seq: &[u8], consensus: &[u8]) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }
    for i in 0..seq.len() {
        let base_mask = base_to_bit(seq[i]);
        let cons_mask = iupac_to_mask(consensus[i]);
        if base_mask & cons_mask == 0 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitmask_roundtrip() {
        let codes = b"ACGTRYSWKMBDHVN";
        for &code in codes {
            let mask = iupac_to_mask(code);
            assert_eq!(
                IUPAC_FROM_MASK[mask as usize], code,
                "Roundtrip failed for '{}'", code as char
            );
        }
    }

    #[test]
    fn test_base_to_bit() {
        assert_eq!(base_to_bit(b'A'), 0b0001);
        assert_eq!(base_to_bit(b'C'), 0b0010);
        assert_eq!(base_to_bit(b'G'), 0b0100);
        assert_eq!(base_to_bit(b'T'), 0b1000);
        assert_eq!(base_to_bit(b'X'), 0);
    }

    #[test]
    fn test_sequence_matches_consensus_bytes() {
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"ACGT"));
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"NCGT"));
        assert!(sequence_matches_consensus_bytes(b"ACGT", b"RCGT"));
        assert!(!sequence_matches_consensus_bytes(b"ACGT", b"YCGT"));
        assert!(!sequence_matches_consensus_bytes(b"ACG", b"ACGT"));
    }
}
