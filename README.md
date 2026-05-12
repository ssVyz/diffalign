# diffalign

A command-line tool for differential oligonucleotide screening. Given a
template DNA sequence and a set of reference sequences, it walks every oligo
window across the template, aligns it against each reference, groups the
matched fragments into distinct variants, and records how many variants are
needed to reach a target coverage threshold. Optionally, it also evaluates
each window against an exclusivity (off-target) FASTA and reports a mismatch
histogram per position.

The output is a single JSON file conforming to [`results_format.md`](results_format.md).

---

## Building

Requires Rust (stable, edition 2021). Build the release binary:

```
cargo build --release
```

The binary lands at `target/release/diffalign` (or `diffalign.exe` on
Windows). Place it wherever you like; on first run it will look for its
config file alongside itself.

---

## Quick start

```
# 1. Create the config file next to the binary.
diffalign --mkini

# 2. Run a screening.
diffalign template.fasta references.fasta -o results.json

# 3. Run with differential / exclusivity mode.
diffalign template.fasta references.fasta -d offtargets.fasta -o results.json
```

File extensions are not enforced. `template.fa`, `template.fasta`,
`template.txt`, or simply `template` are all valid as long as the contents
parse as FASTA. The output path is used verbatim — `-o results` writes a
file literally named `results` (no extension is appended).

By default the tool will **refuse to overwrite an existing output file** —
either delete it, or pick a different path.

---

## Synopsis

```
diffalign <TEMPLATE> <REFERENCE> [-d <EXCLUSIVITY>] -o <OUTPUT> [options]
diffalign --mkini
diffalign --check-config
```

### Required arguments

| Argument             | Description                                                 |
|----------------------|-------------------------------------------------------------|
| `<TEMPLATE>`         | Template FASTA (single sequence; `A`, `C`, `G`, `T` only).  |
| `<REFERENCE>`        | Reference FASTA (multi-sequence).                           |
| `-o, --output PATH`  | Output JSON file. Errors out if the path already exists.    |

### Optional inputs

| Flag                          | Description                                                                 |
|-------------------------------|-----------------------------------------------------------------------------|
| `-d, --diff <FASTA>`          | Exclusivity / off-target FASTA. Presence of this flag enables differential mode. |

### Analysis method

The variant-finding algorithm. Default lives in the INI; CLI flags override
it.

| Flag                                | Default | Description                                                              |
|-------------------------------------|---------|--------------------------------------------------------------------------|
| `--method none\|fixed\|incremental` | `none`  | Pick the variant-finding strategy. See [How it works](#how-it-works).    |
| `--fixed-ambiguities N`             | `1`     | Max IUPAC ambiguity codes per consensus when `--method fixed`.           |
| `--incremental-pct N`               | `50`    | Target coverage % per step when `--method incremental`.                  |
| `--incremental-max-amb N\|none`     | `none`  | Max ambiguities per consensus for `incremental`. `none` = unlimited.     |

### Aligner backend

| Flag                                      | Default    | Description                                                                                  |
|-------------------------------------------|------------|----------------------------------------------------------------------------------------------|
| `-a, --aligner pairwise\|simple\|simple_simd` | `pairwise` | Pick the alignment backend. See [How it works](#how-it-works) for what each does.            |

- `pairwise` — rust-bio Smith-Waterman local alignment. Forward strand only. No oligo-length limit.
- `simple` — bitap (substitutions-only). Forward **and** reverse-complement strands; when the best hit is on the reverse strand the matched fragment is reverse-complemented before being grouped, so downstream variant analysis sees a consistent orientation. Oligo length **must be ≤ 64 bp**.
- `simple_simd` — same algorithm as `simple`, but the inner reference loop is AVX2-vectorized across 4 references per CPU core. Output is bit-identical to `simple` on the same inputs. Requires AVX2: the program detects the CPU at startup and aborts with a message pointing at `--aligner simple` if AVX2 is unavailable. On non-x86_64 builds, the kind is rejected at the same point.

### Pairwise alignment (only used when `--aligner pairwise`)

| Flag                       | Default | Description                                              |
|----------------------------|---------|----------------------------------------------------------|
| `--match-score N`          | `2`     | Score for a matching base.                               |
| `--mismatch-score N`       | `-1`    | Score for a substitution.                                |
| `--gap-open-penalty N`     | `-2`    | Gap-open penalty.                                        |
| `--gap-extend-penalty N`   | `-1`    | Gap-extend penalty.                                      |
| `--max-mismatches N`       | `8`     | Reject a reference alignment above this many mismatches. Also applied to the bitap backends — single flag, shared. |

### Window / range

| Flag                       | Default | Description                                                              |
|----------------------------|---------|--------------------------------------------------------------------------|
| `--min-oligo-length N`     | `18`    | Shortest oligo length to test.                                           |
| `--max-oligo-length N`     | `25`    | Longest oligo length to test.                                            |
| `--length-skip N`          | `0`     | Lengths to skip between processed lengths. `1` = every other length.     |
| `--resolution N`           | `1`     | Position step within a length window. `2` = every other position.        |
| `--coverage-threshold N`   | `90.0`  | Target cumulative variant coverage percentage (0-100).                   |

### Behavior

| Flag                              | Description                                                                   |
|-----------------------------------|-------------------------------------------------------------------------------|
| `--exclude-n` / `--no-exclude-n`  | Include or exclude consensus variants that would require the `N` ambiguity.   |
| `--threads-percent N`             | Percentage of available CPU cores to use (1-100). Floor, min 1 thread.        |
| `-q, --quiet`                     | Suppress the progress bar.                                                    |
| `--timer`                         | Print elapsed wall-clock time at the end of the run.                          |
| `--config PATH`                   | Override the INI file location.                                               |

### Config management

| Flag              | Description                                                                          |
|-------------------|--------------------------------------------------------------------------------------|
| `--mkini`         | Write a default `diffalign.ini` next to the binary and exit. Refuses to overwrite.   |
| `--check-config`  | Report whether a usable INI is present at the expected path. Exit 2 if missing/invalid. |

---

## Configuration file

`diffalign` looks for `diffalign.ini` **next to its own binary**. The file
holds defaults for every parameter; CLI flags override anything set there.

If the file is missing, on an interactive terminal the tool will offer to
create it; otherwise it errors out with a message pointing at `--mkini`.

### Generated default

```ini
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

; Maximum number of variants recorded per position. Leave empty (or 0) for
; no limit. When the limit is reached, the remaining variants' counts are
; folded into the no-match category for that position.
var_limit =

[aligner]
; Which alignment backend to use: pairwise | simple | simple_simd
;   pairwise    = rust-bio Smith-Waterman local alignment
;   simple      = bitap (substitutions-only). Max oligo length is 64 bp;
;                 scans both forward and reverse-complement strands.
;   simple_simd = same algorithm as simple, AVX2-vectorized across references.
;                 Requires a CPU with AVX2; the program errors out at startup
;                 if AVX2 is not detected. Output is bit-identical to simple.
kind = pairwise

[pairwise]
match_score = 2
mismatch_score = -1
gap_open_penalty = -2
gap_extend_penalty = -1
max_mismatches = 8

[simple]
; Maximum substitutions permitted per match for the bitap backends.
max_mismatches = 8

[threads]
; Percentage of available CPU cores to use (1-100). Floor with min of 1 thread.
percent = 100
```

### Precedence

CLI flag > INI file > built-in default.

---

## How it works

### 1. Window enumeration

For each oligo length `L` in `[min_oligo_length, max_oligo_length]` (stepping
by `length_skip + 1`), the tool slides a window of length `L` along the
template. Within a length, positions advance by `resolution`.

### 2. Alignment

For each window, `diffalign` aligns the oligo against every reference
sequence using the configured aligner backend. A reference is counted as a
**match** only if the alignment:

- covers the full oligo (no soft-clipping),
- contains no gaps (no insertions or deletions), and
- has at most `max_mismatches` substitutions.

References that fail any of these checks are tallied as `no_match_count` and
contribute zero coverage.

The matched fragment extracted from each reference (gap-free, length `L`) is
the input to the variant-finding step.

Three backends are available, selected via `--aligner`:

- **`pairwise`** (default) — rust-bio Smith-Waterman local alignment with the
  scoring parameters in the `[pairwise]` INI section. Forward strand only.
  No oligo-length limit. Gapped or partially-covered alignments are computed
  but rejected at the accept layer above.

- **`simple`** — single-`u64` bitap (Wu-Manber substitutions-only variant).
  Both forward and reverse-complement strands are scanned. When the best hit
  is on the reverse strand, the matched fragment is reverse-complemented
  before being passed to variant analysis so all fragments share the oligo's
  orientation. **Oligo length must be ≤ 64 bp.** Lower constant overhead per
  reference than Smith-Waterman.

- **`simple_simd`** — algorithm-identical to `simple`, but the per-reference
  scan loop is AVX2-vectorized across 4 lanes (one reference per lane) on
  x86_64. Output is bit-identical to `simple` on the same inputs (verified by
  a regression test in `screener.rs`). The program checks for AVX2 at startup
  and aborts with a clear message if the CPU does not advertise it; on
  non-x86_64 build targets the kind is rejected up front.

Mismatch-count parity between `pairwise` and the bitap backends is documented
but not exact in edge cases involving IUPAC codes, `N`, or `-` in references:
Smith-Waterman scores those as one mismatch under the configured
`mismatch_score`, while bitap treats any non-ACGT byte as a forced mismatch
at that position. Accept/reject decisions agree in nearly all cases.

### 3. Variant finding

Three methods are available:

- **`none`** (`NoAmbiguities`) — group identical matched fragments and count
  them. The output contains every distinct exact variant, sorted by count
  descending. No IUPAC ambiguity codes are emitted.

- **`fixed`** (`FixedAmbiguities(N)`) — greedy set cover. Repeatedly pick a
  consensus that covers the most still-uncovered sequences while introducing
  at most `N` IUPAC ambiguity codes. Each chosen consensus becomes one
  variant.

- **`incremental`** (`Incremental(P, max_amb)`) — at each step, find a
  consensus that covers at least `P%` of the *currently remaining* sequences,
  optionally bounded by `max_amb` ambiguities. Sequences matched by that
  consensus are removed and the process repeats until none remain.

Variant percentages are computed against `total_sequences` (i.e., including
unmatched references), so unmatched references reduce coverage.

### 4. Coverage threshold

For each window, the tool counts how many top-ranked variants are needed for
the cumulative coverage to reach `coverage_threshold` and records both
`variants_for_threshold` and the actual `coverage_at_threshold` achieved.

### 5. Differential mode (optional)

If `--diff` is given, `diffalign` performs a second pass per window: it
aligns the template oligo against every exclusivity sequence and bins the
results by mismatch count. Each bin records the count and an example name.
The histogram is sorted ascending by mismatch count; sequences that fail to
align (per the same gap/coverage/`max_mismatches` rules above) are bucketed
under the sentinel `mismatches = 4294967295`. Lower mismatch counts indicate
worse exclusivity (the oligo also matches off-target).

### 6. Parallelism

Within each oligo length, positions are processed in parallel via `rayon`,
with one alignment scratch buffer per worker thread. The number of workers
is `floor(threads_percent * available_cores / 100)`, with a minimum of 1.

### 7. Output

A single pretty-printed JSON file. The structure is fully documented in
[`results_format.md`](results_format.md). Highlights:

- `results_by_length` is keyed by oligo length and ordered ascending.
- The `thread_count` recorded in the file is always the concrete count
  actually used for the run (`{ "Fixed": N }`).
- `length_skip` is omitted from output when it is `0`, so default-config
  output stays byte-identical to files produced before this option existed.

---

## Examples

### Default screening, 18-25 bp oligos, no ambiguities

```
diffalign target.fasta refs.fasta -o results.json
```

### Screen 20 bp oligos only, every other position, 50% threads, with off-target check

```
diffalign target.fasta refs.fasta -d offtargets.fasta \
  --min-oligo-length 20 --max-oligo-length 20 \
  --resolution 2 --threads-percent 50 \
  -o results.json
```

### Bitap (faster) with AVX2 vectorization

```
diffalign target.fasta refs.fasta -a simple_simd -o results.json
```

Runs the bitap backend with AVX2 vectorization across references. Requires a
CPU with AVX2; aborts at startup with an error if not available. Output is
bit-identical to `-a simple`. Oligo length must be ≤ 64 bp.

### Length sweep with skip and incremental method

```
diffalign target.fasta refs.fasta \
  --min-oligo-length 18 --max-oligo-length 28 --length-skip 1 \
  --method incremental --incremental-pct 60 --incremental-max-amb 2 \
  --coverage-threshold 95 \
  -o results.json
```

This processes lengths 18, 20, 22, 24, 26, 28; at each window it finds
consensus variants covering at least 60% of remaining sequences (with up to
2 IUPAC codes each) and records how many are needed to reach 95% cumulative
coverage.

---

## Exit codes

| Code | Meaning                                                              |
|------|----------------------------------------------------------------------|
| `0`  | Success.                                                             |
| `1`  | Error (bad input, missing/invalid INI without `--check-config`, etc.).|
| `2`  | `--check-config` reported the INI as missing or invalid.             |

---

## License

See repository root.
