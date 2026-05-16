//! Data types for oligo analysis.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use super::simplescreen::screener::Orientation;

/// Analysis method selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMethod {
    NoAmbiguities,
    FixedAmbiguities(u32),
    /// (target_percentage, optional_max_ambiguities)
    Incremental(u32, Option<u32>),
}

impl Default for AnalysisMethod {
    fn default() -> Self {
        Self::NoAmbiguities
    }
}

impl AnalysisMethod {
    pub fn description(&self) -> String {
        match self {
            Self::NoAmbiguities => "No Ambiguities (exact variants only)".to_string(),
            Self::FixedAmbiguities(n) => format!("Fixed Ambiguities (max {} per variant)", n),
            Self::Incremental(pct, _) => format!("Incremental ({}% coverage per step)", pct),
        }
    }
}

/// Thread count configuration as serialized in the result file.
///
/// The CLI never selects `Auto`; it always resolves the user-supplied
/// percentage to a concrete `Fixed(N)` so the result file records the actual
/// number of cores used. `Auto` is retained only for compatibility with files
/// produced by the original Tauri tool.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThreadCount {
    Auto,
    Fixed(usize),
}

impl Default for ThreadCount {
    fn default() -> Self {
        Self::Auto
    }
}

impl ThreadCount {
    pub fn get_count(&self) -> usize {
        match self {
            Self::Auto => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            Self::Fixed(n) => *n,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PairwiseParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
    pub max_mismatches: u32,
}

impl Default for PairwiseParams {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open_penalty: -2,
            gap_extend_penalty: -1,
            max_mismatches: 8,
        }
    }
}

/// Parameters for the simplescreen (bitap) aligner.
///
/// `max_mismatches` is the only knob the bitap algorithm needs. The four
/// scoring fields on `PairwiseParams` are SW-specific and don't apply here.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SimpleParams {
    pub max_mismatches: u32,
}

impl Default for SimpleParams {
    fn default() -> Self {
        Self { max_mismatches: 8 }
    }
}

/// Which alignment backend to use for screening.
///
/// `snake_case` rename is identical to `lowercase` for the single-word variants
/// (`Pairwise`, `Simple`) so legacy JSON files continue to round-trip; it only
/// affects `SimpleSimd` → `"simple_simd"` and `SimpleCuda` → `"simple_cuda"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlignerKind {
    Pairwise,
    Simple,
    SimpleSimd,
    /// GPU bitap. Only present in builds compiled with the `cuda` feature;
    /// kept as a serde-recognized variant on all builds so JSON files that
    /// mention `simple_cuda` still parse (the CLI rejects it at startup on
    /// non-cuda builds).
    SimpleCuda,
}

impl Default for AlignerKind {
    fn default() -> Self {
        Self::Pairwise
    }
}

impl AlignerKind {
    pub fn is_default(&self) -> bool {
        matches!(self, Self::Pairwise)
    }

    /// True if this kind uses the bitap (simplescreen) algorithm under the
    /// hood (whether scalar, SIMD-vectorized, or GPU-accelerated).
    pub fn is_bitap(&self) -> bool {
        matches!(self, Self::Simple | Self::SimpleSimd | Self::SimpleCuda)
    }
}

fn is_zero_u32(n: &u32) -> bool {
    *n == 0
}

fn is_false(b: &bool) -> bool {
    !*b
}

fn is_default_aligner(k: &AlignerKind) -> bool {
    k.is_default()
}

/// Per-reference anchor position produced by a single alignment pass.
///
/// Used internally by anchored screening to record where the anchor oligo
/// matched in each reference, so length-`L` fragments can be derived from
/// that same position later instead of re-aligning per length. Not serialized.
#[derive(Debug, Clone, Copy)]
pub struct AnchorHit {
    /// 0-based start position on the forward strand of the reference.
    pub start: usize,
    pub orientation: Orientation,
    pub mismatches: u32,
}

pub const DEFAULT_MAX_SEEDS: u32 = 50;

fn default_max_seeds() -> u32 {
    DEFAULT_MAX_SEEDS
}

fn is_default_max_seeds(n: &u32) -> bool {
    *n == DEFAULT_MAX_SEEDS
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisParams {
    pub method: AnalysisMethod,
    pub pairwise: PairwiseParams,
    /// Selected alignment backend. Omitted from JSON output when set to the
    /// default (`Pairwise`) so default-config output stays byte-identical to
    /// files written before the simple-aligner option existed.
    #[serde(default, skip_serializing_if = "is_default_aligner")]
    pub aligner: AlignerKind,
    /// Parameters for the simple (bitap) aligner. Only present in JSON output
    /// when `aligner = Simple`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub simple: Option<SimpleParams>,
    pub exclude_n: bool,
    pub min_oligo_length: u32,
    pub max_oligo_length: u32,
    pub resolution: u32,
    pub coverage_threshold: f64,
    pub thread_count: ThreadCount,
    /// Number of lengths to skip between processed lengths within
    /// `[min_oligo_length, max_oligo_length]`. `0` (default) processes every
    /// length; `1` processes every other length, and so on.
    /// Omitted from JSON output when zero so default-config output stays
    /// byte-identical to files written by the original tool.
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub length_skip: u32,
    /// Maximum number of variants recorded per position. `None` = unlimited.
    /// When set and exceeded, the trailing variants are dropped and their
    /// counts are folded into the no-match category for that position.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub var_limit: Option<u32>,
    /// Number of seed sequences the ambiguity-aware variant finder tries per
    /// greedy step. Higher = better coverage at linear cost; default 50. Only
    /// used by `FixedAmbiguities` and `Incremental` methods.
    #[serde(default = "default_max_seeds", skip_serializing_if = "is_default_max_seeds")]
    pub max_seeds: u32,
    /// When true, run the per-position search once at `anchored_length` and
    /// reuse those positions to derive matched fragments for every length in
    /// `[min_oligo_length, max_oligo_length]`. Defaults to false (legacy
    /// behavior). Omitted from JSON output when false.
    #[serde(default, skip_serializing_if = "is_false")]
    pub anchored: bool,
    /// Length of the anchor search when `anchored` is true. `None` falls back
    /// to `min_oligo_length`. Must lie within `[min_oligo_length,
    /// max_oligo_length]`. Omitted from JSON output when `None`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchored_length: Option<u32>,
}

impl Default for AnalysisParams {
    fn default() -> Self {
        Self {
            method: AnalysisMethod::NoAmbiguities,
            pairwise: PairwiseParams::default(),
            aligner: AlignerKind::Pairwise,
            simple: None,
            exclude_n: true,
            min_oligo_length: 18,
            max_oligo_length: 25,
            resolution: 1,
            coverage_threshold: 90.0,
            thread_count: ThreadCount::Auto,
            length_skip: 0,
            var_limit: None,
            max_seeds: DEFAULT_MAX_SEEDS,
            anchored: false,
            anchored_length: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variant {
    pub sequence: String,
    pub count: usize,
    pub percentage: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowAnalysisResult {
    pub variants: Vec<Variant>,
    pub total_sequences: usize,
    pub sequences_analyzed: usize,
    pub no_match_count: usize,
    pub variants_for_threshold: usize,
    pub coverage_at_threshold: f64,
    pub skipped: bool,
    pub skip_reason: Option<String>,
}

impl Default for WindowAnalysisResult {
    fn default() -> Self {
        Self {
            variants: Vec::new(),
            total_sequences: 0,
            sequences_analyzed: 0,
            no_match_count: 0,
            variants_for_threshold: 0,
            coverage_at_threshold: 0.0,
            skipped: false,
            skip_reason: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LengthResult {
    pub oligo_length: u32,
    pub positions: Vec<PositionResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionResult {
    pub position: usize,
    pub variants_needed: usize,
    pub analysis: WindowAnalysisResult,
    #[serde(default)]
    pub exclusivity: Option<ExclusivityResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusivityResult {
    pub total_sequences: usize,
    pub no_match_count: usize,
    pub mismatch_histogram: Vec<MismatchBucket>,
    pub min_mismatches: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MismatchBucket {
    pub mismatches: u32,
    pub count: usize,
    pub example_name: String,
}

/// User annotation on the template sequence.
///
/// Carried through for output compatibility with files produced by the
/// original Tauri program. The CLI always emits an empty annotations array.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Annotation {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub direction: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreeningResults {
    pub params: AnalysisParams,
    pub template_length: usize,
    pub total_sequences: usize,
    pub template_sequence: String,
    /// `BTreeMap` keyed by oligo length so JSON key order is stable across runs.
    pub results_by_length: BTreeMap<u32, LengthResult>,
    #[serde(default)]
    pub differential_enabled: bool,
    #[serde(default)]
    pub exclusivity_sequence_count: Option<usize>,
    #[serde(default)]
    pub annotations: Vec<Annotation>,
}

impl ScreeningResults {
    pub fn new(
        params: AnalysisParams,
        template_length: usize,
        total_sequences: usize,
        template_sequence: String,
        differential_enabled: bool,
        exclusivity_sequence_count: Option<usize>,
    ) -> Self {
        Self {
            params,
            template_length,
            total_sequences,
            template_sequence,
            results_by_length: BTreeMap::new(),
            differential_enabled,
            exclusivity_sequence_count,
            annotations: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ProgressUpdate {
    pub current_length: u32,
    pub current_position: usize,
    pub total_positions: usize,
    pub lengths_completed: u32,
    pub total_lengths: u32,
    pub message: String,
}
