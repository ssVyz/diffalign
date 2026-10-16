# diffalign

A command-line tool for differential oligonucleotide screening. Given a
template DNA sequence and a set of reference sequences, it walks every oligo
window across the template, aligns it against each reference, groups the
matched fragments into distinct variants, and records how many variants are
needed to reach a target coverage threshold. Optionally, it also evaluates
each window against an exclusivity (off-target) FASTA and reports a mismatch
histogram per position.

The output is a single JSON file conforming to [`results_format.md`](results_format.md).


## QUICKSTART GUIDE

*if you have no time and you just want to run it for testing*

**Step 1:** Clone diffalign from github if you haven't already. (run from a folder of your choice)

> git clone https://github.com/ssVyz/diffalign

... then: get the result viewer

> git clone https://github.com/ssVyz/diffalign_result_viewer

... then: get the example sequence set

> git clone https://github.com/ssVyz/diffalign_example_set


**Step 2:** Compile the binary. You need a Rust toolchain. You also need Python for the result viewer if you are using it.
... build from the project folder:

> cargo build --release

... binary will appear in your chosen folder ../diffalign/target/release
... either add to path or copy/paste to wherever you want to run it.

**Step 3:** Run diffalign on the example files.

> diffalign flua_pdm09_m_refseq.fasta FluA_m_1y_0806_sequences.fasta -o test_output.json

... at the location of the test files. It will emit a test_output.json as your result

**Step 4:** Open the file with the result viewer. You need astral UV, or manually install libraries.
... in the result viewer root directory:

> uv sync

> uv run main.py

... Then load your results via the GUI.


---

## AI disclosure

This tool and its precursors were created via agentic coding using various versions
of Claude code (Anthropic) and chatGPT (OpenAI). Large parts of this readme file were
written by the coding agents.

## Building

Requires Rust (stable, edition 2021). Build the release binary:

```
cargo build --release
```

The binary lands at `target/release/diffalign` (or `diffalign.exe` on
Windows). Place it wherever you like. It looks for an optional config file
alongside itself and uses built-in defaults when there is none (see
[Configuration file](#configuration-file)).

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

## Use (simplified)

```
# 1. Run a screening with the built-in defaults.
diffalign template.fasta references.fasta -o results.json

# 2. Run with differential / exclusivity mode.
diffalign template.fasta references.fasta -d offtargets.fasta -o results.json

# 3. Optional: write the defaults to diffalign.ini next to the binary, to edit them.
diffalign --mkini
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

The variant-finding algorithm. The **Default** columns in this section show
the built-in defaults, which are also what `--mkini` writes. Values in
`diffalign.ini` override them, and CLI flags override both.

| Flag                                | Default | Description                                                              |
|-------------------------------------|---------|--------------------------------------------------------------------------|
| `--method none\|fixed\|incremental` | `none`  | Pick the variant-finding strategy. See [How it works](#how-it-works).    |
| `--fixed-ambiguities N`             | `1`     | Max IUPAC ambiguity codes per consensus when `--method fixed`.           |
| `--incremental-pct N`               | `66`    | Target coverage % per step when `--method incremental` (1-100).          |
| `--incremental-max-amb N\|none`     | `none`  | Max ambiguities per consensus for `incremental`. `none` = unlimited.     |
| `--max-seeds N`                     | `1000`  | Seed sequences tried per greedy step by `fixed` and `incremental`. Higher = better coverage, linearly slower. Must be >= 1. |
| `--var-limit N\|none`               | `15`    | Maximum variants recorded per position. Variants beyond the limit are dropped and their references counted as no-match. `none` (or `0`) = unlimited. |

### Aligner backend

| Flag                                                         | Default    | Description                                                                                  |
|--------------------------------------------------------------|------------|----------------------------------------------------------------------------------------------|
| `-a, --aligner pairwise\|simple\|simple_simd\|simple_cuda`   | `simple`   | Pick the alignment backend. See [How it works](#how-it-works) for what each does.            |
| `--max-mismatches N`                                         | `10` (bitap) / `8` (pairwise) | Maximum substitutions in an accepted hit. The INI keeps separate values in `[simple]` and `[pairwise]`; this flag sets both. With a bitap backend it must be smaller than `--min-oligo-length`. |

All backends search both the forward and the reverse-complement strand. When
the best hit is on the reverse strand, the matched fragment is
reverse-complemented before being grouped, so variant analysis sees a
consistent orientation.

- `pairwise` — rust-bio Smith-Waterman local alignment. No oligo-length limit.
- `simple` — bitap (substitutions-only). Oligo length **must be ≤ 64 bp**.
- `simple_simd` — same algorithm as `simple`, but the inner reference loop is AVX2-vectorized across 4 references per CPU core. Output is bit-identical to `simple` on the same inputs. Requires AVX2: the program detects the CPU at startup and aborts with a message pointing at `--aligner simple` if AVX2 is unavailable. On non-x86_64 builds, the kind is rejected at the same point.
- `simple_cuda` — same algorithm as `simple`, but the per-reference scan runs on an NVIDIA GPU with one CUDA thread per reference. Output is bit-identical to `simple` on the same inputs, and works with every `--method` (the GPU only does the alignment stage; variant analysis runs on the CPU afterwards). Only available in builds compiled with `--features cuda`. Caps `max_mismatches` at 16. Requires a CUDA-capable GPU and the CUDA runtime at startup; the program errors out cleanly if either is missing.

### Pairwise alignment (only used when `--aligner pairwise`)

| Flag                       | Default | Description                                              |
|----------------------------|---------|----------------------------------------------------------|
| `--match-score N`          | `2`     | Score for a matching base.                               |
| `--mismatch-score N`       | `-1`    | Score for a substitution.                                |
| `--gap-open-penalty N`     | `-2`    | Gap-open penalty.                                        |
| `--gap-extend-penalty N`   | `-1`    | Gap-extend penalty.                                      |

The mismatch cap for `pairwise` is set with `--max-mismatches` (see
[Aligner backend](#aligner-backend)).

### Window / range

| Flag                       | Default | Description                                                              |
|----------------------------|---------|--------------------------------------------------------------------------|
| `--min-oligo-length N`     | `16`    | Shortest oligo length to test.                                           |
| `--max-oligo-length N`     | `30`    | Longest oligo length to test (≤ 64 for the bitap backends).              |
| `--length-skip N`          | `0`     | Lengths to skip between processed lengths. `1` = every other length.     |
| `--resolution N`           | `1`     | Position step within a length window. `2` = every other position.        |
| `--coverage-threshold N`   | `95.0`  | Target cumulative variant coverage percentage (0-100).                   |

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
| `--exclude-n` / `--no-exclude-n`  | Forbid (default) or allow the `N` code in consensus variants built by `fixed` and `incremental`. Does not remove reference sequences. |
| `--threads-percent N`             | Percentage of available CPU cores to use (1-100, default 100). Floor, min 1 thread. |
| `-q, --quiet`                     | Suppress the progress bar.                                                    |
| `--timer`                         | Print elapsed wall-clock time at the end of the run.                          |
| `--config PATH`                   | Use this INI file instead of the one next to the binary. The file must exist. |

### Config management

| Flag              | Description                                                                          |
|-------------------|--------------------------------------------------------------------------------------|
| `--mkini`         | Write a default `diffalign.ini` next to the binary and exit. Refuses to overwrite.   |
| `--check-config`  | Report whether a usable INI is present at the expected path. Exit 2 if missing/invalid (the program can run without an ini file, which then use the built-in defaults. Outside of testing, a correctly configured ini file will be required for meaningful results). |

---

## Configuration file

`diffalign` looks for `diffalign.ini` **next to its own binary**. The file
is optional: it lets you change the defaults for every parameter, and CLI
flags override anything set there.

If the file is missing, `diffalign` uses its built-in defaults and prints a
one-line note. These are the same values `diffalign --mkini` writes (shown
below), so running without an INI behaves exactly like running with a
freshly generated one. Keys left out of an INI also fall back to the
built-in defaults. A file named explicitly with `--config` must exist.

### Generated default

This is the file `diffalign --mkini` writes:

```ini
; diffalign.ini — default configuration
;
; CLI flags always override values in this file.

[analysis]
; Variant-finding method: none | fixed | incremental
method = none

; Used when method = fixed: maximum ambiguity codes per consensus variant.
fixed_ambiguities = 1

; Used when method = incremental: target coverage percentage per step (1-100).
incremental_pct = 66

; Used when method = incremental: maximum ambiguity codes per variant.
; Leave empty for unlimited.
incremental_max_ambiguities =

; Whether to exclude variants that would require an N ambiguity code.
exclude_n = true

; Oligo length range (inclusive).
min_oligo_length = 16
max_oligo_length = 30

; Number of lengths to skip between processed lengths.
; 0 = process every length, 1 = every other (e.g. 20, 22, 24), etc.
length_skip = 0

; Position step within a length window. 1 = every position.
resolution = 1

; Target cumulative variant coverage percentage (0-100).
; Required for backwards compatibility.
coverage_threshold = 95.0

; Maximum number of variants recorded per position.
; Leave empty (or 0) for no limit. When the limit is reached, the
; remaining variants' counts are folded into the no-match category
; for that position. Useful for keeping output JSON files manageable.
var_limit = 15

; Number of seed sequences the ambiguity-aware variant finder tries
; per greedy step. Only used by methods `fixed` and `incremental`.
; Higher values explore more starting points (better coverage) at
; linear runtime cost. Must be >= 1.
max_seeds = 1000

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
kind = simple

[pairwise]
match_score = 2
mismatch_score = -1
gap_open_penalty = -2
gap_extend_penalty = -1
max_mismatches = 8

[simple]
; Maximum substitutions permitted per match for the bitap backend.
max_mismatches = 10

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
sequence, on both strands, using the configured aligner backend. A reference
is counted as a **match** only if its best alignment:

- covers the full oligo (no soft-clipping),
- contains no gaps (no insertions or deletions), and
- has at most `max_mismatches` substitutions (the `[simple]` value for the
  bitap backends, the `[pairwise]` value for `pairwise`).

References that fail any of these checks are tallied as `no_match_count` and
contribute zero coverage.

The matched fragment extracted from each reference (gap-free, length `L`) is
the input to the variant-finding step.

Four backends are available, selected via `--aligner`:

- **`pairwise`** — rust-bio Smith-Waterman local alignment with the
  scoring parameters in the `[pairwise]` INI section. The oligo and its
  reverse complement are both aligned; the accepted hit with fewer
  mismatches wins, and a reverse-strand fragment is reverse-complemented.
  No oligo-length limit. Gapped or partially-covered alignments are computed
  but rejected at the accept layer above.

- **`simple`** (default) — single-`u64` bitap (Wu-Manber substitutions-only variant).
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
  descending (ties alphabetically). This method adds no IUPAC ambiguity
  codes; any that appear were copied verbatim from a reference sequence.

- **`fixed`** (`FixedAmbiguities(N)`) — greedy set cover. Repeatedly pick a
  consensus that covers the most still-uncovered sequences while introducing
  at most `N` IUPAC ambiguity codes. Each chosen consensus becomes one
  variant, listed in the order it was picked.

- **`incremental`** (`Incremental(P, max_amb)`) — at each step, find the
  consensus with the fewest ambiguity codes that covers at least `P%` of the
  *currently remaining* sequences, optionally bounded by `max_amb`
  ambiguities. If none reaches `P%`, the best-covering consensus found is
  used. Sequences matched by that consensus are removed and the process
  repeats until none remain. Variants are listed in the order they were
  found, so a later variant can have a higher count than an earlier one.

`fixed` and `incremental` try up to `max_seeds` seed sequences per step and,
with `exclude_n` (the default), never use the `N` code. Every matched
reference is counted in exactly one variant: the first one that covers it.

Variant percentages are computed against `total_sequences` (i.e., including
unmatched references), so unmatched references reduce coverage.

If `var_limit` is set (default 15), only the first `var_limit` variants of a
window are kept. The references belonging to dropped variants are counted as
no-match.

### 4. Coverage threshold

For each window, the tool walks the variant list in the order described
above and counts how many variants are needed for the cumulative coverage to
reach `coverage_threshold`. It records both `variants_for_threshold` and the
actual `coverage_at_threshold` achieved. If the threshold cannot be reached
(too many unmatched references, or variants dropped by `var_limit`), all
listed variants are counted and `coverage_at_threshold` stays below the
threshold.

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

Because the anchor only runs at template positions where the anchor-length
window fits, lengths shorter than the anchor stop at position
`template_length - anchored_length` instead of `template_length - L`.

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
- Some parameters (e.g. `length_skip` when it is `0`, or `aligner` when it is
  `pairwise`) are left out of the output when they hold the value older
  versions used implicitly, so older result files keep parsing. See
  [`results_format.md`](results_format.md) for the full list.

---

## Examples

### Default screening

```
diffalign target.fasta refs.fasta -o results.json
```

Uses the built-in defaults (or your `diffalign.ini`): 16-30 bp oligos, bitap
aligner with up to 10 mismatches, no ambiguities, 95% coverage threshold,
at most 15 variants per position.

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
--features cuda`. Works with every `--method`; aborts at startup with an
error if a usable GPU + CUDA runtime is not available. Output is
bit-identical to `-a simple`. Oligo length must be ≤ 64 bp;
`max_mismatches` must be ≤ 16 (the default of 10 is fine). Differential
mode (`-d offtargets.fasta`) works the same way as with the other aligners.

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
| `1`  | Error (bad input, invalid INI, missing `--config` file, etc.).       |
| `2`  | `--check-config` reported the INI as missing or invalid.             |

---

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
