# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.4] - 2026-05-29

### Fixed
- `--method` on the command line now correctly inherits the INI's per-method
  parameters for the newly selected method. Previously, switching the method
  via CLI (e.g. `--method incremental` against an INI whose `method = none`)
  silently discarded the INI's `incremental_pct` / `incremental_max_ambiguities`
  (and equivalently for `fixed_ambiguities`) and fell back to built-in defaults
  — 50% target coverage and *unlimited* ambiguity codes per variant. In
  practice this produced "incremental" variants with far more ambiguities than
  the configured cap (predominantly Ns when `exclude_n = false`, other IUPAC
  codes when `exclude_n = true`). Each method's parameters are now stored on
  the loaded config independently of which `method` is active, so a CLI
  override picks them up.

### Added
- Regression tests for the CLI/INI method merge:
  `config::tests::retains_incremental_keys_when_method_is_none` pins that
  per-method INI keys are preserved when another method is active, and
  `cli::tests::cli_method_override_inherits_ini_incremental_params` pins the
  end-to-end behavior of `--method incremental` on top of `method = none`.

## [0.9.3] - 2026-05-21

### Fixed
- Pairwise aligner (`--aligner pairwise`) now scans both the forward strand
  and the reverse complement of each reference, matching the long-standing
  behavior of the bitap aligners (`simple` / `simple_simd` / `simple_cuda`).
  Previously, references where the target oligo appeared only on the reverse
  strand were silently counted as `no_match` under the pairwise backend. This
  affected `variants` / `count` / `percentage` / `variants_for_threshold` in
  the main analysis, the `mismatch_histogram` in differential mode (RC-only
  off-targets were misclassified as `u32::MAX` "no match"), and anchor
  derivation in anchored mode. When the reverse-strand pass wins, the
  matched fragment is reverse-complemented back to oligo orientation before
  variant grouping, and `AnchorHit.orientation` is recorded correctly.

### Added
- Strand-handling parity tests:
  `pairwise::reverse_complement_match_returns_oligo_oriented_fragment`,
  `pairwise::forward_preferred_over_reverse_on_ties`,
  `pairwise::anchor_orientation_recorded_for_reverse_hits`, and the
  end-to-end `screener::pairwise_finds_reverse_complement_matches` so
  the divergence cannot regress unnoticed.

## [0.9.2] - 2026-05-21

### Changed
- Rewrote `results_format.md` to match the JSON actually produced by the
  current `diffalign` binary. Replaced legacy "Oligoscreen Desktop" branding,
  documented the new `params` fields (`aligner`, `simple`, `length_skip`,
  `var_limit`, `max_seeds`, `anchored`, `anchored_length`), documented the
  `annotations` array on `ScreeningResults`, expanded `AlignerKind` to
  include `simple_simd` and `simple_cuda`, clarified `ThreadCount` (the CLI
  always emits `Fixed`; `Auto` is read-only compatibility), and added
  examples for bitap + differential mode and for anchored mode.

## [0.9.1] - 2026-05-21

### Added
- `CHANGELOG.md` and de-ignored it (along with `CLAUDE.md`) in `.gitignore`.
- `CLAUDE.md` with versioning workflow instructions for future Claude sessions.

### Changed
- Bumped `Cargo.toml` version from `0.1.0` to `0.9.1` to reflect the road-to-1.0
  state described in `version1_todo.md`.
