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

### Optional: CUDA backend

The GPU aligner (`--aligner simple_cuda`) is gated behind a Cargo feature.
To build with it enabled:

```
cargo build --release --features cuda
```

Build requirements (on top of Rust):

- NVIDIA CUDA Toolkit (tested with 13.2; minimum supported GPU is Turing /
  sm_75 on CUDA 13.x). The toolkit is found via `CUDA_PATH`, falling back to
  the standard install location on Windows.
- A C++ host compiler that `nvcc` recognizes — on Windows that's MSVC
  (Visual Studio Build Tools or full Visual Studio), located automatically
  via the `cc` crate.

Runtime requirements when actually using `--aligner simple_cuda`:

- An NVIDIA GPU of supported compute capability (sm_75 / Turing or newer).
- The CUDA runtime library (`cudart64_*.dll` on Windows, `libcudart.so` on
  Linux) on the system search path.

The compiled binary still works on systems without a GPU — only
`--aligner simple_cuda` requires the runtime. The other aligners run
normally and the program errors out cleanly with a helpful message if the
GPU backend is selected but unusable.

A plain `cargo build --release` (no `--features cuda`) produces a smaller
binary with no CUDA runtime dependency at all.

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

| Flag                                                         | Default    | Description                                                                                  |
|--------------------------------------------------------------|------------|----------------------------------------------------------------------------------------------|
| `-a, --aligner pairwise\|simple\|simple_simd\|simple_cuda`   | `pairwise` | Pick the alignment backend. See [How it works](#how-it-works) for what each does.            |

- `pairwise` — rust-bio Smith-Waterman local alignment. Forward strand only. No oligo-length limit.
- `simple` — bitap (substitutions-only). Forward **and** reverse-complement strands; when the best hit is on the reverse strand the matched fragment is reverse-complemented before being grouped, so downstream variant analysis sees a consistent orientation. Oligo length **must be ≤ 64 bp**.
- `simple_simd` — same algorithm as `simple`, but the inner reference loop is AVX2-vectorized across 4 references per CPU core. Output is bit-identical to `simple` on the same inputs. Requires AVX2: the program detects the CPU at startup and aborts with a message pointing at `--aligner simple` if AVX2 is unavailable. On non-x86_64 builds, the kind is rejected at the same point.
- `simple_cuda` — same algorithm as `simple`, but the per-reference scan runs on an NVIDIA GPU with one CUDA thread per reference. Output is bit-identical to `simple` on the same inputs, and works with every `--method` (the GPU only does the alignment stage; variant analysis runs on the CPU afterwards). Only available in builds compiled with `--features cuda`. Caps `max_mismatches` at 16. Requires a CUDA-capable GPU and the CUDA runtime at startup; the program errors out cleanly if either is missing.

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

### Anchored mode

Opt-in mode that runs the per-position alignment once at a single anchor
length and derives every other length's matched fragments from those stored
positions, instead of re-running the search per length. Much faster on wide
length ranges; the position bias of the anchor length carries through to
every derived length.

| Flag                  | Default                 | Description                                                       |
|-----------------------|-------------------------|-------------------------------------------------------------------|
| `--anchored`          | off                     | Enable anchored mode.                                             |
| `--anchored-length N` | `--min-oligo-length`    | Length used for the anchor search. Must lie within `[min, max]`.  |

See [Anchored mode](#anchored-mode-details) below for the semantics.

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

; Anchored mode: when true, the per-position search runs once at
; `anchored_length`, and every length in the length range derives its
; variants from those stored positions (instead of re-running the search
; per length).
anchored = false

; Length used for the anchor search when anchored = true. Must lie within
; [min_oligo_length, max_oligo_length]. Leave empty to fall back to
; min_oligo_length.
anchored_length =

[aligner]
; Which alignment backend to use: pairwise | simple | simple_simd | simple_cuda
;   pairwise    = rust-bio Smith-Waterman local alignment
;   simple      = bitap (substitutions-only). Max oligo length is 64 bp;
;                 scans both forward and reverse-complement strands.
;   simple_simd = same algorithm as simple, AVX2-vectorized across references.
;                 Requires a CPU with AVX2; the program errors out at startup
;                 if AVX2 is not detected. Output is bit-identical to simple.
;   simple_cuda = same algorithm as simple, GPU-accelerated. Only available
;                 in builds compiled with --features cuda. Works with every
;                 method. Requires an NVIDIA GPU + CUDA runtime;
;                 caps max_mismatches at 16. Output is bit-identical to simple.
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

- **`simple_cuda`** — algorithm-identical to `simple`, runs on an NVIDIA GPU
  with one CUDA thread per reference. References are uploaded to the GPU
  once at run start; per-window kernel launches reuse the same on-device
  buffers. Output is bit-identical to `simple` on the same inputs (regression
  test in `screener.rs`, gated on `cuda` feature + a usable GPU). Only
  available in builds compiled with `--features cuda` (see
  [Building](#building)). The GPU performs only the alignment stage, so every
  variant-finding method works as it does for `simple`; `max_mismatches` is
  capped at 16 (a kernel register-pressure limit). The program checks for a
  usable CUDA device at startup and refuses to run with a clear error if the
  runtime, driver, or GPU is missing; the other aligners are unaffected on
  systems without CUDA.

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

### 5b. Anchored mode (optional, opt-in) <a id="anchored-mode-details"></a>

By default, the per-position search runs once per oligo length, which is
unbiased but scales linearly with the number of lengths in
`[min_oligo_length, max_oligo_length]`. With `--anchored`, the search runs
**once at `--anchored-length`** (default = `--min-oligo-length`), and the
matched per-reference position from that single pass is reused to derive the
length-`L` matched fragments for every length in the range.

For an anchor matched at reference position `s` with anchor length `L_a`,
the length-`L` fragment is extracted as follows:

- **Forward orientation:** `reference[s..s+L]` — same 5'-start as the
  anchor; right-extend for `L > L_a`, right-truncate for `L < L_a`.
- **Reverse orientation:** `reverse_complement(reference[s+L_a-L..s+L_a])`
  — the same convention, applied to the reverse strand and RC-flipped back
  into oligo orientation.

After extraction, the length-`L` fragment is re-checked against the template
oligo `template[p..p+L]` using a base-by-base Hamming distance
(case-insensitive ACGT; any other byte counts as a mismatch, matching the
bitap convention). If the Hamming distance exceeds `max_mismatches`, the
reference is no-match **for that length only** — the anchor itself, and
other lengths that still fit, are unaffected. The same is true if the
extension/truncation would run past the end of the reference.

If no anchor was found for a given `(template_position, reference)` pair,
the reference is no-match at **every** length for that position.

The same logic governs the exclusivity (differential) pass under `-d`: the
anchor is run once on the exclusivity set, and the mismatch histogram for
each length is computed from the Hamming distances at the anchor's stored
positions.

Constraints: `anchored_length` must lie within
`[min_oligo_length, max_oligo_length]`. Without an explicit value it
defaults to `min_oligo_length`. The mode works with every `--aligner`
backend and every `--method`.

Trade-off: anchored mode trades flexibility for speed. The default mode
re-searches per length and can find a different best position per length;
anchored mode locks each reference's position to whatever the anchor
length found. Use the default when bias matters; use anchored when you
need a fast wide-length sweep.

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

### Bitap on GPU (largest reference sets)

```
diffalign target.fasta refs.fasta -a simple_cuda -o results.json
```

Runs the bitap backend on an NVIDIA GPU with one CUDA thread per reference.
Only available when the binary was built with `cargo build --release
--features cuda`. Restricted to `--method none`; aborts at startup with an
error if a usable GPU + CUDA runtime is not available, or if combined with
an unsupported method. Output is bit-identical to `-a simple`. Oligo length
must be ≤ 64 bp; `max_mismatches` must be ≤ 16. Differential mode
(`-d offtargets.fasta`) works the same way as with the other aligners.

### Anchored mode: fast wide-length sweep

```
diffalign target.fasta refs.fasta -a simple_simd \
  --min-oligo-length 18 --max-oligo-length 30 \
  --anchored --anchored-length 18 \
  -o results.json
```

Runs the per-position search once at length 18 and derives the matched
fragments for lengths 19..30 from those stored positions, instead of
re-running the search per length. `--anchored-length` defaults to
`--min-oligo-length` and may be omitted (`--anchored` alone is enough).
References whose extension would run past the reference end, or whose
length-`L` fragment exceeds `max_mismatches` against the template, become
no-match for *that length only*. Works with every aligner and method.

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
