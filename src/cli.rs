//! Command-line argument parsing and merge with the INI config.

use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::{ArgAction, Parser};

use crate::analysis::{
    AlignerKind, AnalysisMethod, AnalysisParams, PairwiseParams, SimpleParams, ThreadCount,
};
use crate::config::{Config, parse_aligner_kind};

#[derive(Debug, Parser)]
#[command(
    name = "diffalign",
    about = "Differential oligonucleotide screening (CLI)",
    version,
    long_about = "Screens a template DNA sequence against reference sequences \
to identify primer/probe sites, with optional differential screening \
against an exclusivity (off-target) FASTA. CLI flags override values \
from diffalign.ini."
)]
pub struct Cli {
    /// Template FASTA (single sequence; A/C/G/T only).
    #[arg(required_unless_present_any = ["mkini", "check_config"])]
    pub template: Option<PathBuf>,

    /// Reference FASTA (multi-sequence).
    #[arg(required_unless_present_any = ["mkini", "check_config"])]
    pub reference: Option<PathBuf>,

    /// Differential / exclusivity FASTA. Triggers differential mode.
    #[arg(short = 'd', long = "diff", value_name = "FASTA")]
    pub diff: Option<PathBuf>,

    /// Output JSON file. Errors out if the file already exists.
    #[arg(
        short = 'o',
        long = "output",
        value_name = "PATH",
        required_unless_present_any = ["mkini", "check_config"]
    )]
    pub output: Option<PathBuf>,

    // ── analysis method ───────────────────────────────────────────────
    /// Variant-finding method: none | fixed | incremental.
    #[arg(long, value_parser = ["none", "fixed", "incremental"])]
    pub method: Option<String>,

    /// Max ambiguity codes per consensus (used with --method fixed).
    #[arg(long)]
    pub fixed_ambiguities: Option<u32>,

    /// Target coverage % per step (used with --method incremental).
    #[arg(long)]
    pub incremental_pct: Option<u32>,

    /// Max ambiguity codes per variant (used with --method incremental).
    /// Pass `none` for unlimited.
    #[arg(long, value_name = "N|none")]
    pub incremental_max_amb: Option<String>,

    // ── aligner backend ───────────────────────────────────────────────
    /// Alignment backend: pairwise (Smith-Waterman) | simple (bitap, ≤64 bp) |
    /// simple_simd (AVX2-vectorized bitap; requires CPU AVX2 support).
    #[arg(short = 'a', long = "aligner", value_parser = ["pairwise", "simple", "simple_simd", "simd"])]
    pub aligner: Option<String>,

    // ── pairwise ──────────────────────────────────────────────────────
    #[arg(long)]
    pub match_score: Option<i32>,
    #[arg(long)]
    pub mismatch_score: Option<i32>,
    #[arg(long)]
    pub gap_open_penalty: Option<i32>,
    #[arg(long)]
    pub gap_extend_penalty: Option<i32>,
    /// Maximum substitutions allowed per match. Applies to whichever
    /// aligner backend is selected.
    #[arg(long)]
    pub max_mismatches: Option<u32>,

    // ── exclude_n ─────────────────────────────────────────────────────
    /// Exclude variants that would require an N ambiguity.
    #[arg(long, action = ArgAction::SetTrue, conflicts_with = "no_exclude_n")]
    pub exclude_n: bool,
    /// Allow N ambiguity in variants.
    #[arg(long = "no-exclude-n", action = ArgAction::SetTrue)]
    pub no_exclude_n: bool,

    // ── lengths / windows ─────────────────────────────────────────────
    #[arg(long)]
    pub min_oligo_length: Option<u32>,
    #[arg(long)]
    pub max_oligo_length: Option<u32>,
    /// Number of lengths to skip between processed lengths.
    /// 0 = every length, 1 = every other (e.g. 20, 22, 24).
    #[arg(long)]
    pub length_skip: Option<u32>,
    /// Position step within a length window. 1 = every position.
    #[arg(long)]
    pub resolution: Option<u32>,
    #[arg(long)]
    pub coverage_threshold: Option<f64>,

    /// Maximum number of variants recorded per position.
    /// Pass `none` (or `0`) for no limit. When set and exceeded, the
    /// remaining variants' counts are folded into the no-match category
    /// for that position.
    #[arg(long = "var-limit", value_name = "N|none")]
    pub var_limit: Option<String>,

    // ── threads ───────────────────────────────────────────────────────
    /// Percentage of available CPU cores to use (1-100).
    /// Resolved to a concrete core count (floor, min 1) at run time.
    #[arg(long)]
    pub threads_percent: Option<u32>,

    // ── progress ──────────────────────────────────────────────────────
    /// Suppress the progress bar.
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Print elapsed wall-clock time at the end of the run.
    #[arg(long)]
    pub timer: bool,

    // ── config management ─────────────────────────────────────────────
    /// Write a default diffalign.ini next to the binary and exit.
    #[arg(long, exclusive = true)]
    pub mkini: bool,

    /// Print whether a usable diffalign.ini is present and exit.
    #[arg(long, exclusive = true)]
    pub check_config: bool,

    /// Override the location of the INI file (advanced).
    #[arg(long, value_name = "PATH")]
    pub config: Option<PathBuf>,
}

impl Cli {
    fn exclude_n_override(&self) -> Option<bool> {
        if self.exclude_n {
            Some(true)
        } else if self.no_exclude_n {
            Some(false)
        } else {
            None
        }
    }

    /// Apply CLI overrides on top of `cfg` and produce the final
    /// `AnalysisParams` that get fed into the screening engine.
    ///
    /// Resolves the threads percentage to a concrete `ThreadCount::Fixed(N)`
    /// based on the number of cores available right now.
    pub fn into_analysis_params(&self, cfg: &Config) -> Result<AnalysisParams> {
        // Build method, allowing per-flag overrides on top of the INI choice.
        let method = self.resolve_method(cfg)?;

        let aligner = match &self.aligner {
            Some(s) => parse_aligner_kind(s)?,
            None => cfg.aligner,
        };

        let pairwise = PairwiseParams {
            match_score: self.match_score.unwrap_or(cfg.pairwise.match_score),
            mismatch_score: self.mismatch_score.unwrap_or(cfg.pairwise.mismatch_score),
            gap_open_penalty: self.gap_open_penalty.unwrap_or(cfg.pairwise.gap_open_penalty),
            gap_extend_penalty: self.gap_extend_penalty.unwrap_or(cfg.pairwise.gap_extend_penalty),
            max_mismatches: self.max_mismatches.unwrap_or(cfg.pairwise.max_mismatches),
        };

        let simple = SimpleParams {
            max_mismatches: self.max_mismatches.unwrap_or(cfg.simple.max_mismatches),
        };

        let exclude_n = self.exclude_n_override().unwrap_or(cfg.exclude_n);

        let min_oligo_length = self.min_oligo_length.unwrap_or(cfg.min_oligo_length);
        let max_oligo_length = self.max_oligo_length.unwrap_or(cfg.max_oligo_length);
        let length_skip = self.length_skip.unwrap_or(cfg.length_skip);
        let resolution = self.resolution.unwrap_or(cfg.resolution);
        let coverage_threshold = self.coverage_threshold.unwrap_or(cfg.coverage_threshold);
        let var_limit = match &self.var_limit {
            Some(s) => parse_optional_u32(s, "--var-limit")?.filter(|&n| n > 0),
            None => cfg.var_limit,
        };

        if min_oligo_length == 0 {
            bail!("min_oligo_length must be >= 1");
        }
        if max_oligo_length < min_oligo_length {
            bail!(
                "max_oligo_length ({}) must be >= min_oligo_length ({})",
                max_oligo_length,
                min_oligo_length
            );
        }
        if resolution == 0 {
            bail!("resolution must be >= 1");
        }
        if !(0.0..=100.0).contains(&coverage_threshold) {
            bail!(
                "coverage_threshold must be between 0 and 100 (got {})",
                coverage_threshold
            );
        }
        if aligner.is_bitap() {
            const SIMPLE_MAX: u32 = 64;
            if max_oligo_length > SIMPLE_MAX {
                bail!(
                    "aligner = {} supports oligos up to {} bp; got max_oligo_length = {}.\n\
                     Lower max_oligo_length (and min_oligo_length) to <= {}, or use --aligner pairwise.",
                    aligner_name(aligner),
                    SIMPLE_MAX,
                    max_oligo_length,
                    SIMPLE_MAX
                );
            }
            if simple.max_mismatches >= min_oligo_length {
                bail!(
                    "max_mismatches ({}) must be strictly less than min_oligo_length ({}) when aligner = {}",
                    simple.max_mismatches,
                    min_oligo_length,
                    aligner_name(aligner)
                );
            }
        }

        if aligner == AlignerKind::SimpleSimd {
            ensure_avx2_available()?;
        }

        let threads_percent = self.threads_percent.unwrap_or(cfg.threads_percent);
        if threads_percent == 0 || threads_percent > 100 {
            bail!(
                "threads percent must be between 1 and 100 (got {})",
                threads_percent
            );
        }
        let thread_count = resolve_threads(threads_percent);

        let simple_field = if aligner.is_bitap() {
            Some(simple)
        } else {
            None
        };

        Ok(AnalysisParams {
            method,
            pairwise,
            aligner,
            simple: simple_field,
            exclude_n,
            min_oligo_length,
            max_oligo_length,
            resolution,
            coverage_threshold,
            thread_count,
            length_skip,
            var_limit,
        })
    }

    fn resolve_method(&self, cfg: &Config) -> Result<AnalysisMethod> {
        // Defaults from the INI...
        let (mut method_name, mut fixed_amb, mut incr_pct, mut incr_max_amb): (
            &str,
            u32,
            u32,
            Option<u32>,
        ) = match cfg.method {
            AnalysisMethod::NoAmbiguities => ("none", 1, 50, None),
            AnalysisMethod::FixedAmbiguities(n) => ("fixed", n, 50, None),
            AnalysisMethod::Incremental(p, m) => ("incremental", 1, p, m),
        };

        // ...overridden per-field by CLI flags.
        if let Some(ref m) = self.method {
            method_name = m.as_str();
        }
        if let Some(v) = self.fixed_ambiguities {
            fixed_amb = v;
        }
        if let Some(v) = self.incremental_pct {
            incr_pct = v;
        }
        if let Some(v) = &self.incremental_max_amb {
            incr_max_amb = parse_optional_u32(v, "--incremental-max-amb")?;
        }

        match method_name {
            "none" | "no" | "no_ambiguities" => Ok(AnalysisMethod::NoAmbiguities),
            "fixed" => Ok(AnalysisMethod::FixedAmbiguities(fixed_amb)),
            "incremental" => {
                if !(1..=100).contains(&incr_pct) {
                    bail!("incremental percentage must be between 1 and 100");
                }
                Ok(AnalysisMethod::Incremental(incr_pct, incr_max_amb))
            }
            other => bail!("unknown method '{}'", other),
        }
    }
}

fn aligner_name(k: AlignerKind) -> &'static str {
    match k {
        AlignerKind::Pairwise => "pairwise",
        AlignerKind::Simple => "simple",
        AlignerKind::SimpleSimd => "simple_simd",
    }
}

/// Refuse to run when `simple_simd` was selected but the CPU (or build target)
/// does not advertise AVX2.
fn ensure_avx2_available() -> Result<()> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            return Ok(());
        }
        bail!(
            "aligner = simple_simd requires AVX2, but this CPU does not advertise it.\n\
             Use --aligner simple (same algorithm, scalar) or --aligner pairwise instead."
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        bail!(
            "aligner = simple_simd is only available on x86_64 builds.\n\
             Use --aligner simple (same algorithm, scalar) or --aligner pairwise instead."
        );
    }
}

fn parse_optional_u32(value: &str, flag_name: &str) -> Result<Option<u32>> {
    let trimmed = value.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("none") {
        return Ok(None);
    }
    trimmed
        .parse::<u32>()
        .map(Some)
        .map_err(|e| anyhow::anyhow!("invalid value for {}: {}", flag_name, e))
}

/// Resolve a percentage to a concrete `ThreadCount::Fixed(N)`.
/// Floor with a minimum of 1 thread.
pub fn resolve_threads(percent: u32) -> ThreadCount {
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let raw = (cores as u64 * percent as u64) / 100;
    let n = raw.max(1) as usize;
    ThreadCount::Fixed(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_threads_floor_min_1() {
        // Whatever the host has, 1% should still resolve to >= 1.
        let tc = resolve_threads(1);
        match tc {
            ThreadCount::Fixed(n) => assert!(n >= 1),
            _ => panic!("expected Fixed"),
        }
    }

    #[test]
    fn resolve_threads_100_percent_uses_all_cores() {
        let cores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let tc = resolve_threads(100);
        match tc {
            ThreadCount::Fixed(n) => assert_eq!(n, cores),
            _ => panic!("expected Fixed"),
        }
    }

    #[test]
    fn parse_optional_u32_handles_none_keyword() {
        assert_eq!(parse_optional_u32("none", "--x").unwrap(), None);
        assert_eq!(parse_optional_u32("NONE", "--x").unwrap(), None);
        assert_eq!(parse_optional_u32("", "--x").unwrap(), None);
        assert_eq!(parse_optional_u32("3", "--x").unwrap(), Some(3));
        assert!(parse_optional_u32("foo", "--x").is_err());
    }
}
