//! Analysis engine for oligonucleotide screening.
//!
//! Lifted from the reference Tauri program; only data-structure and parameter
//! changes were made. Behavior of the screening, alignment, and variant logic
//! is unchanged.

#![allow(dead_code)]

mod analyzer;
#[cfg(feature = "cuda")]
mod cuda_align;
mod fasta;
mod iupac;
mod pairwise;
mod screener;
mod simd_align;
mod simple_align;
mod simplescreen;
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
#[allow(unused_imports)]
pub use simple_align::*;
#[allow(unused_imports)]
pub use simd_align::*;
#[cfg(feature = "cuda")]
#[allow(unused_imports)]
pub use cuda_align::*;
