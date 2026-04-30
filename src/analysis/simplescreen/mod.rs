//! Vendored simplescreen sources (bitap, pattern, iupac, screener).
//!
//! The four files in this folder are kept as close to upstream as possible —
//! the only edits are `crate::` → `super::` so they compile as siblings under
//! `analysis::simplescreen`. See `simplescreen-src/simplescreen_README.md` in
//! the upstream tool for the original layout. Updating the algorithm is a
//! copy-paste of the four files plus those four import patches.

#![allow(dead_code)]

pub mod bitap;
pub mod iupac;
pub mod pattern;
pub mod screener;
