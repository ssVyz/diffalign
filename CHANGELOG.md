# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
