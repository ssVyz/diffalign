//! Analysis engine for oligonucleotide screening.
//!
//! Lifted from the reference Tauri program; only data-structure and parameter
//! changes were made. Behavior of the screening, alignment, and variant logic
//! is unchanged.

#![allow(dead_code)]

mod analyzer;
mod fasta;
mod iupac;
mod pairwise;
mod screener;
mod types;

pub use fasta::*;
pub use screener::*;
pub use types::*;

#[allow(unused_imports)]
pub use analyzer::*;
#[allow(unused_imports)]
pub use iupac::*;
#[allow(unused_imports)]
pub use pairwise::*;
