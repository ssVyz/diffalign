//! INI configuration file: defaults for analysis parameters.
//!
//! Lives next to the binary as `diffalign.ini`. CLI flags override values
//! loaded from this file. A default file can be written via `--mkini`.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};
use ini::Ini;

use crate::analysis::{AlignerKind, AnalysisMethod, DEFAULT_MAX_SEEDS, PairwiseParams, SimpleParams};

pub const INI_FILE_NAME: &str = "diffalign.ini";

/// Fully-resolved configuration loaded from an INI file (or built-in defaults).
///
/// Mirrors `AnalysisParams` but stores threads as a percentage of available
/// cores; the percentage is resolved to a concrete core count at run time.
#[derive(Debug, Clone)]
pub struct Config {
    pub method: AnalysisMethod,
    pub aligner: AlignerKind,
    pub pairwise: PairwiseParams,
    pub simple: SimpleParams,
    pub exclude_n: bool,
    pub min_oligo_length: u32,
    pub max_oligo_length: u32,
    pub length_skip: u32,
    pub resolution: u32,
    pub coverage_threshold: f64,
    /// Threads as a percentage of available cores. 0 < x <= 100.
    pub threads_percent: u32,
    /// Maximum variants recorded per position. `None` = unlimited.
    pub var_limit: Option<u32>,
    /// Number of seed sequences the ambiguity-aware variant finder tries per
    /// greedy step. Only used by `fixed` and `incremental` methods.
    pub max_seeds: u32,
    /// Enable anchored mode: run the per-position search once at
    /// `anchored_length` and derive every length's variants from those
    /// positions. Defaults to false.
    pub anchored: bool,
    /// Length of the anchor search when `anchored = true`. `None` falls back
    /// to `min_oligo_length` at run time.
    pub anchored_length: Option<u32>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            method: AnalysisMethod::NoAmbiguities,
            aligner: AlignerKind::Pairwise,
            pairwise: PairwiseParams::default(),
            simple: SimpleParams::default(),
            exclude_n: true,
            min_oligo_length: 18,
            max_oligo_length: 25,
            length_skip: 0,
            resolution: 1,
            coverage_threshold: 90.0,
            threads_percent: 100,
            var_limit: None,
            max_seeds: DEFAULT_MAX_SEEDS,
            anchored: false,
            anchored_length: None,
        }
    }
}

/// Default INI text written by `--mkini`.
pub const DEFAULT_INI_TEMPLATE: &str = r#"; diffalign.ini — default configuration
;
; Values here are the defaults used when no matching CLI flag is given.
; CLI flags always override values in this file.

[analysis]
; Variant-finding method: none | fixed | incremental
method = none

; Used when method = fixed: maximum ambiguity codes per consensus variant.
fixed_ambiguities = 1

; Used when method = incremental: target coverage percentage per step (1-100).
incremental_pct = 50

; Used when method = incremental: maximum ambiguity codes per variant.
; Leave empty for unlimited.
incremental_max_ambiguities =

; Whether to exclude variants that would require an N ambiguity code.
exclude_n = true

; Oligo length range (inclusive).
min_oligo_length = 18
max_oligo_length = 25

; Number of lengths to skip between processed lengths.
; 0 = process every length, 1 = every other (e.g. 20, 22, 24), etc.
length_skip = 0

; Position step within a length window. 1 = every position.
resolution = 1

; Target cumulative variant coverage percentage (0-100).
coverage_threshold = 90.0

; Maximum number of variants recorded per position.
; Leave empty (or 0) for no limit. When the limit is reached, the
; remaining variants' counts are folded into the no-match category
; for that position. Useful for keeping output JSON files manageable.
var_limit =

; Number of seed sequences the ambiguity-aware variant finder tries
; per greedy step. Only used by methods `fixed` and `incremental`.
; Higher values explore more starting points (better coverage) at
; linear runtime cost; 50 is a good default. Must be >= 1.
max_seeds = 50

; Anchored mode: when true, the per-position search runs once at
; `anchored_length`, and every length in the length range derives its
; variants from those stored positions (instead of re-running the search
; per length). Much faster on wide length ranges; can introduce slight
; bias because the position is fixed by the anchor length. Default false.
anchored = false

; Length used for the anchor search when anchored = true. Must lie within
; [min_oligo_length, max_oligo_length]. Leave empty to fall back to
; min_oligo_length.
anchored_length =

[aligner]
; Which alignment backend to use: pairwise | simple | simple_simd | simple_cuda
;   pairwise    = rust-bio Smith-Waterman local alignment (gaps allowed by score,
;                 gapped/partial-coverage alignments are still rejected at the
;                 accept layer)
;   simple      = bitap (substitutions-only). Faster; max oligo length is 64 bp;
;                 scans both forward and reverse-complement strands.
;   simple_simd = same algorithm as simple, AVX2-vectorized across references.
;                 Requires a CPU with AVX2; the program errors out at startup
;                 if AVX2 is not detected. Output is bit-identical to simple.
;   simple_cuda = same algorithm as simple, GPU-accelerated across references.
;                 Only available in builds compiled with the `cuda` feature.
;                 Works with every method; caps max_mismatches at 16. Requires
;                 an NVIDIA GPU + CUDA runtime at startup; the program errors
;                 out otherwise. Output is bit-identical to simple.
kind = pairwise

[pairwise]
match_score = 2
mismatch_score = -1
gap_open_penalty = -2
gap_extend_penalty = -1
max_mismatches = 8

[simple]
; Maximum substitutions permitted per match for the bitap backend.
max_mismatches = 8

[threads]
; Percentage of available CPU cores to use (1-100). Floor with min of 1 thread.
percent = 100
"#;

/// Resolve the on-disk location of the INI file (next to the binary).
pub fn ini_path() -> Result<PathBuf> {
    let exe = std::env::current_exe().context("failed to locate the diffalign executable")?;
    let dir = exe
        .parent()
        .ok_or_else(|| anyhow!("executable path has no parent directory"))?;
    Ok(dir.join(INI_FILE_NAME))
}

/// Write the default INI to `path`. Errors if the file already exists.
pub fn write_default_ini(path: &Path) -> Result<()> {
    if path.exists() {
        bail!("INI file already exists at {}", path.display());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("creating directory {}", parent.display()))?;
        }
    }
    fs::write(path, DEFAULT_INI_TEMPLATE)
        .with_context(|| format!("writing default INI to {}", path.display()))?;
    Ok(())
}

/// Load and validate a `Config` from an INI file.
pub fn load(path: &Path) -> Result<Config> {
    if !path.exists() {
        bail!("INI file not found at {}", path.display());
    }

    let ini = Ini::load_from_file(path)
        .with_context(|| format!("parsing INI file {}", path.display()))?;

    let mut cfg = Config::default();

    // [analysis]
    let analysis = ini.section(Some("analysis"));
    let method_str = get_str(analysis, "method").unwrap_or_else(|| "none".to_string());
    let fixed_amb = get_u32(analysis, "fixed_ambiguities")?.unwrap_or(1);
    let incr_pct = get_u32(analysis, "incremental_pct")?.unwrap_or(50);
    let incr_max_amb = get_optional_u32(analysis, "incremental_max_ambiguities")?;
    cfg.method = parse_method(&method_str, fixed_amb, incr_pct, incr_max_amb)?;
    cfg.exclude_n = get_bool(analysis, "exclude_n")?.unwrap_or(cfg.exclude_n);
    cfg.min_oligo_length = get_u32(analysis, "min_oligo_length")?.unwrap_or(cfg.min_oligo_length);
    cfg.max_oligo_length = get_u32(analysis, "max_oligo_length")?.unwrap_or(cfg.max_oligo_length);
    cfg.length_skip = get_u32(analysis, "length_skip")?.unwrap_or(cfg.length_skip);
    cfg.resolution = get_u32(analysis, "resolution")?.unwrap_or(cfg.resolution);
    cfg.coverage_threshold =
        get_f64(analysis, "coverage_threshold")?.unwrap_or(cfg.coverage_threshold);
    // var_limit: empty or 0 means unlimited.
    cfg.var_limit = get_optional_u32(analysis, "var_limit")?.filter(|&n| n > 0);
    cfg.max_seeds = get_u32(analysis, "max_seeds")?.unwrap_or(cfg.max_seeds);
    cfg.anchored = get_bool(analysis, "anchored")?.unwrap_or(cfg.anchored);
    cfg.anchored_length = get_optional_u32(analysis, "anchored_length")?;

    // [aligner]
    let aligner = ini.section(Some("aligner"));
    if let Some(kind) = get_str(aligner, "kind") {
        cfg.aligner = parse_aligner_kind(&kind)?;
    }

    // [pairwise]
    let pairwise = ini.section(Some("pairwise"));
    cfg.pairwise.match_score =
        get_i32(pairwise, "match_score")?.unwrap_or(cfg.pairwise.match_score);
    cfg.pairwise.mismatch_score =
        get_i32(pairwise, "mismatch_score")?.unwrap_or(cfg.pairwise.mismatch_score);
    cfg.pairwise.gap_open_penalty =
        get_i32(pairwise, "gap_open_penalty")?.unwrap_or(cfg.pairwise.gap_open_penalty);
    cfg.pairwise.gap_extend_penalty =
        get_i32(pairwise, "gap_extend_penalty")?.unwrap_or(cfg.pairwise.gap_extend_penalty);
    cfg.pairwise.max_mismatches =
        get_u32(pairwise, "max_mismatches")?.unwrap_or(cfg.pairwise.max_mismatches);

    // [simple]
    let simple = ini.section(Some("simple"));
    cfg.simple.max_mismatches =
        get_u32(simple, "max_mismatches")?.unwrap_or(cfg.simple.max_mismatches);

    // [threads]
    let threads = ini.section(Some("threads"));
    cfg.threads_percent = get_u32(threads, "percent")?.unwrap_or(cfg.threads_percent);

    cfg.validate()?;
    Ok(cfg)
}

pub fn parse_aligner_kind(name: &str) -> Result<AlignerKind> {
    match name.trim().to_ascii_lowercase().as_str() {
        "pairwise" | "pw" | "rustbio" => Ok(AlignerKind::Pairwise),
        "simple" | "simplescreen" | "bitap" => Ok(AlignerKind::Simple),
        "simple_simd" | "simd" | "simplesimd" => Ok(AlignerKind::SimpleSimd),
        "simple_cuda" | "cuda" | "gpu" | "simplecuda" => Ok(AlignerKind::SimpleCuda),
        other => bail!(
            "unknown aligner kind '{}' (expected: pairwise | simple | simple_simd | simple_cuda)",
            other
        ),
    }
}

impl Config {
    pub fn validate(&self) -> Result<()> {
        if self.min_oligo_length == 0 {
            bail!("min_oligo_length must be >= 1");
        }
        if self.max_oligo_length < self.min_oligo_length {
            bail!(
                "max_oligo_length ({}) must be >= min_oligo_length ({})",
                self.max_oligo_length,
                self.min_oligo_length
            );
        }
        if self.resolution == 0 {
            bail!("resolution must be >= 1");
        }
        if !(0.0..=100.0).contains(&self.coverage_threshold) {
            bail!(
                "coverage_threshold must be between 0 and 100 (got {})",
                self.coverage_threshold
            );
        }
        if self.threads_percent == 0 || self.threads_percent > 100 {
            bail!(
                "threads percent must be between 1 and 100 (got {})",
                self.threads_percent
            );
        }
        if self.max_seeds == 0 {
            bail!("max_seeds must be >= 1");
        }
        Ok(())
    }
}

fn parse_method(
    name: &str,
    fixed_amb: u32,
    incr_pct: u32,
    incr_max_amb: Option<u32>,
) -> Result<AnalysisMethod> {
    match name.trim().to_ascii_lowercase().as_str() {
        "none" | "no" | "no_ambiguities" => Ok(AnalysisMethod::NoAmbiguities),
        "fixed" => Ok(AnalysisMethod::FixedAmbiguities(fixed_amb)),
        "incremental" => {
            if !(1..=100).contains(&incr_pct) {
                bail!(
                    "incremental_pct must be between 1 and 100 (got {})",
                    incr_pct
                );
            }
            Ok(AnalysisMethod::Incremental(incr_pct, incr_max_amb))
        }
        other => bail!(
            "unknown analysis method '{}' (expected: none | fixed | incremental)",
            other
        ),
    }
}

fn get_str(section: Option<&ini::Properties>, key: &str) -> Option<String> {
    section
        .and_then(|s| s.get(key))
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
}

fn get_u32(section: Option<&ini::Properties>, key: &str) -> Result<Option<u32>> {
    match get_str(section, key) {
        None => Ok(None),
        Some(v) => v
            .parse::<u32>()
            .map(Some)
            .map_err(|e| anyhow!("invalid integer for '{}': {}", key, e)),
    }
}

fn get_optional_u32(section: Option<&ini::Properties>, key: &str) -> Result<Option<u32>> {
    // Same as get_u32 but the empty-string treatment in get_str already maps
    // empty to None, which is exactly what we want for "unlimited".
    get_u32(section, key)
}

fn get_i32(section: Option<&ini::Properties>, key: &str) -> Result<Option<i32>> {
    match get_str(section, key) {
        None => Ok(None),
        Some(v) => v
            .parse::<i32>()
            .map(Some)
            .map_err(|e| anyhow!("invalid integer for '{}': {}", key, e)),
    }
}

fn get_f64(section: Option<&ini::Properties>, key: &str) -> Result<Option<f64>> {
    match get_str(section, key) {
        None => Ok(None),
        Some(v) => v
            .parse::<f64>()
            .map(Some)
            .map_err(|e| anyhow!("invalid number for '{}': {}", key, e)),
    }
}

fn get_bool(section: Option<&ini::Properties>, key: &str) -> Result<Option<bool>> {
    match get_str(section, key) {
        None => Ok(None),
        Some(v) => match v.to_ascii_lowercase().as_str() {
            "true" | "yes" | "1" | "on" => Ok(Some(true)),
            "false" | "no" | "0" | "off" => Ok(Some(false)),
            other => bail!("invalid boolean for '{}': '{}'", key, other),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_default_template() {
        let dir = tempdir();
        let path = dir.join("diffalign.ini");
        fs::write(&path, DEFAULT_INI_TEMPLATE).unwrap();
        let cfg = load(&path).unwrap();
        assert_eq!(cfg.method, AnalysisMethod::NoAmbiguities);
        assert_eq!(cfg.min_oligo_length, 18);
        assert_eq!(cfg.max_oligo_length, 25);
        assert_eq!(cfg.length_skip, 0);
        assert_eq!(cfg.threads_percent, 100);
        assert!(cfg.exclude_n);
    }

    #[test]
    fn rejects_bad_threads_percent() {
        let dir = tempdir();
        let path = dir.join("diffalign.ini");
        fs::write(&path, "[threads]\npercent = 0\n").unwrap();
        assert!(load(&path).is_err());
    }

    #[test]
    fn parses_incremental_method() {
        let dir = tempdir();
        let path = dir.join("diffalign.ini");
        let body = r#"
[analysis]
method = incremental
incremental_pct = 25
incremental_max_ambiguities = 2
"#;
        fs::write(&path, body).unwrap();
        let cfg = load(&path).unwrap();
        assert_eq!(cfg.method, AnalysisMethod::Incremental(25, Some(2)));
    }

    #[test]
    fn parses_incremental_unlimited_amb() {
        let dir = tempdir();
        let path = dir.join("diffalign.ini");
        let body = r#"
[analysis]
method = incremental
incremental_pct = 50
incremental_max_ambiguities =
"#;
        fs::write(&path, body).unwrap();
        let cfg = load(&path).unwrap();
        assert_eq!(cfg.method, AnalysisMethod::Incremental(50, None));
    }

    #[test]
    fn write_default_refuses_to_overwrite() {
        let dir = tempdir();
        let path = dir.join("diffalign.ini");
        fs::write(&path, "junk").unwrap();
        assert!(write_default_ini(&path).is_err());
    }

    fn tempdir() -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "diffalign-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }
}
