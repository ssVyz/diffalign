//! diffalign — differential oligonucleotide screening.
//!
//! The analysis engine is exposed as a library so it can be exercised from
//! tests and reused. The CLI binary (`src/main.rs`) is a thin wrapper around
//! this crate plus the `cli`, `config`, and `progress` modules.

pub mod analysis;
pub mod cli;
pub mod config;
pub mod key_listener;
pub mod pause;
pub mod progress;
