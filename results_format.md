# diffalign Results File Format Specification

This document is the definitive reference for the JSON result file format
produced by `diffalign`. A file that conforms to this specification can be
round-tripped through `serde_json` against the types defined in
`src/analysis/types.rs`. (`diffalign` itself only writes these files; the
types also derive `Deserialize` so other tools can read them.)

## Overview

- **Format**: JSON
- **Encoding**: UTF-8
- **Serialization**: Rust `serde_json::to_string_pretty` (pretty-printed,
  2-space indent)
- **File extension**: `.json` (convention; not enforced)
- **Output path**: chosen by the user via `--output`; the binary has no
  default filename, refuses to overwrite an existing file, and creates
  missing parent directories.
- **Versioning**: the file does not record which `diffalign` version wrote
  it.
- **Size**: grows with (number of lengths) × (number of positions) ×
  (variants per position). `var_limit` caps the last factor.

## Conventions

- **Coordinates** are 0-based. The *template oligo* for oligo length `L` at
  position `p` is `template_sequence[p .. p + L]` (end exclusive).
- **Orientation**: every reference is searched on both strands. When the best
  hit is on the reverse strand, the matched fragment is reverse-complemented
  before analysis. Every variant sequence is therefore in the same
  orientation as the template and has length `L`, so it can be compared
  base-by-base with the template oligo.
- **Percentages** are numbers in `[0, 100]`, not fractions.

## Defaults and omitted fields

Two different notions of "default" apply to `AnalysisParams`:

1. **Configuration defaults** are what a run uses when the user passes no
   overriding CLI flags: the values written by `diffalign --mkini`. The same
   values are built in and used when no `diffalign.ini` exists. They are
   listed in the **`--mkini` default** column of the `AnalysisParams` table.
2. **Serialization omissions**: a few fields are left out of the JSON when
   they hold one specific value, and that value is assumed when they are
   absent on read. This keeps files compatible with older builds and with the
   original Tauri program, which did not have these fields. The omission
   values are fixed in code and are *independent* of the configuration
   defaults. With the `--mkini` defaults, for example, `aligner`, `simple`,
   `var_limit` and `max_seeds` are all written.

| Field | Omitted from output when | Value assumed when absent |
|---|---|---|
| `aligner` | `"pairwise"` | `"pairwise"` |
| `simple` | the aligner is `"pairwise"` | not set (no bitap parameters) |
| `length_skip` | `0` | `0` |
| `var_limit` | unlimited | unlimited |
| `max_seeds` | exactly `50` | `50` |
| `anchored` | `false` | `false` |
| `anchored_length` | not set | not set (the anchor length is `min_oligo_length`) |

All other `AnalysisParams` fields are always written and are **required** on
read (they have no serde default). Outside `AnalysisParams`, the fields
`differential_enabled`, `exclusivity_sequence_count`, `annotations` and
`PositionResult.exclusivity` are always written. On read they may be absent
and then default to `false`, `null`, `[]` and `null` respectively.

## How each window is computed

For each oligo length `L` and position `p`:

1. Take the template oligo `template_sequence[p .. p + L]`.
2. Search every reference sequence for the oligo on both strands with the
   selected aligner (see [`PairwiseParams`](#pairwiseparams) and
   [`SimpleParams`](#simpleparams) for the acceptance rules). Each reference
   either yields one matched fragment (length `L`, in template orientation)
   or is a no-match.
3. Group the matched fragments into variants according to `method`.
4. Compute each variant's percentage of *all* references, apply `var_limit`,
   and determine how many variants are needed to reach `coverage_threshold`.
5. In differential mode, search every exclusivity sequence for the
   **template oligo** (not the variants) the same way, and build a
   histogram of mismatch counts.

**Anchored mode** (`anchored: true`) replaces steps 2 and 5 with a single
search at the anchor length `A` (`anchored_length`, or `min_oligo_length`
when that is not set). For each length `L`, the fragment of length `L`
starting at the anchor hit is extracted on the hit's strand: it is extended
when `L > A` and truncated when `L < A`. Its Hamming distance to the template
oligo is then checked against the mismatch cap. A fragment that runs past the
end of the sequence, or exceeds the cap, counts as a no-match for that
length.

## Top-Level Object: `ScreeningResults`

```jsonc
{
  "params":                      AnalysisParams,
  "template_length":             integer,       // >= 1
  "total_sequences":             integer,       // >= 1
  "template_sequence":           string,        // uppercase DNA: A, C, G, T only
  "results_by_length":           { "<length>": LengthResult, ... },
  "differential_enabled":        boolean,
  "exclusivity_sequence_count":  integer | null,
  "annotations":                 [Annotation, ...]
}
```

| Field | Type | Description |
|---|---|---|
| `params` | `AnalysisParams` | The parameters that were used to generate these results. |
| `template_length` | integer | Length of `template_sequence` in base pairs. |
| `total_sequences` | integer | Number of sequences in the reference FASTA. Every reference is screened; none are filtered out. |
| `template_sequence` | string | The full template DNA sequence. Uppercase characters `A`, `C`, `G`, `T` only. |
| `results_by_length` | object | Map from oligo length (as a **string** key, e.g. `"18"`) to its `LengthResult`. Contains exactly the processed lengths: `min_oligo_length`, then every `length_skip + 1` up to `max_oligo_length` (`max_oligo_length` itself only if it falls on that step). Keys are emitted in ascending numeric order (backed by a `BTreeMap`). |
| `differential_enabled` | boolean | `true` if an exclusivity FASTA was given (`--diff`) and differential analysis was performed. |
| `exclusivity_sequence_count` | integer or `null` | Number of exclusivity sequences used. `null` when differential mode is disabled. |
| `annotations` | array of `Annotation` | User annotations on the template. The diffalign CLI always writes an empty array; the field exists for compatibility with files produced by the original Tauri program. |

> **Key serialization note**: `results_by_length` keys are stringified
> integers (e.g. `"18"`, `"25"`), not bare integer keys. This is because JSON
> object keys must be strings.

---

## `AnalysisParams`

Fields are written in this order:

```jsonc
{
  "method":             AnalysisMethod,
  "pairwise":           PairwiseParams,
  "aligner":            AlignerKind,         // omitted when "pairwise"
  "simple":             SimpleParams,        // present only for bitap aligners
  "exclude_n":          boolean,
  "min_oligo_length":   integer,             // >= 1
  "max_oligo_length":   integer,             // >= min_oligo_length
  "resolution":         integer,             // >= 1
  "coverage_threshold": number,              // 0.0 - 100.0
  "thread_count":       ThreadCount,
  "length_skip":        integer,             // omitted when 0
  "var_limit":          integer,             // omitted when unlimited
  "max_seeds":          integer,             // omitted when exactly 50
  "anchored":           boolean,             // omitted when false
  "anchored_length":    integer              // omitted when not set
}
```

| Field | Type | `--mkini` default | Description |
|---|---|---|---|
| `method` | `AnalysisMethod` | `"NoAmbiguities"` | How matched fragments are grouped into variants. |
| `pairwise` | `PairwiseParams` | see below | Smith-Waterman parameters. Always written, even when a bitap aligner was used (its values were then not used). |
| `aligner` | `AlignerKind` | `"simple"` | Alignment backend used to search references and exclusivity sequences. |
| `simple` | `SimpleParams` | `{"max_mismatches": 10}` | Bitap parameters. Present exactly when `aligner` is `"simple"`, `"simple_simd"` or `"simple_cuda"`. |
| `exclude_n` | boolean | `true` | When `true`, consensus variants built by `FixedAmbiguities` / `Incremental` never use the `N` code: fragments that would need all four bases at one position are not merged. It does **not** remove reference sequences, and it has no effect on `NoAmbiguities`. |
| `min_oligo_length` | integer | `16` | Shortest oligo length tested. |
| `max_oligo_length` | integer | `30` | Longest oligo length tested. At most `64` for bitap aligners. |
| `resolution` | integer | `1` | Step between analyzed positions: positions `0, r, 2r, …`. |
| `coverage_threshold` | number | `95.0` | Target cumulative coverage, as a percentage of all references, used for `variants_for_threshold`. |
| `thread_count` | `ThreadCount` | `{"Fixed": <all cores>}` | Threads used for the run. See [`ThreadCount`](#threadcount). |
| `length_skip` | integer | `0` | Number of lengths skipped between processed lengths. `0` processes every length; `1` every other one (e.g. 18, 20, 22, 24 for 18–25). |
| `var_limit` | integer | `15` | Maximum number of variants kept per position. Absent = unlimited. Extra variants (the tail of the `variants` array) are dropped; their counts move from `sequences_analyzed` to `no_match_count`. |
| `max_seeds` | integer | `1000` | Number of seed sequences tried per greedy step. Only used by `FixedAmbiguities` and `Incremental`; recorded for every method. |
| `anchored` | boolean | `false` | Whether anchored mode was used. See [How each window is computed](#how-each-window-is-computed). |
| `anchored_length` | integer | not set | Anchor length in anchored mode. Only present when `anchored` is `true`. Absent = `min_oligo_length`. Lies within `[min_oligo_length, max_oligo_length]`. |

### `AnalysisMethod`

An externally tagged enum. Exactly one of the following forms:

| Variant | JSON representation | Description |
|---|---|---|
| No ambiguities | `"NoAmbiguities"` | Every distinct matched fragment is its own variant. |
| Fixed ambiguities | `{"FixedAmbiguities": <integer>}` | Greedy cover. Each step picks the consensus with at most N IUPAC ambiguity codes that covers the most still-unassigned fragments, trying up to `max_seeds` seed sequences. N >= 0. The `--mkini` setting is `1`. |
| Incremental | `{"Incremental": [<integer>, <integer or null>]}` | Each step looks for the consensus with the fewest ambiguity codes that covers at least the given percentage (1–100) of the still-unassigned fragments. If none reaches it, the best-covering consensus found is used. The second element caps the ambiguity codes per variant (`null` = unlimited). The `--mkini` settings give `[66, null]`. |

### `PairwiseParams`

```jsonc
{
  "match_score":        integer,   // --mkini default: 2
  "mismatch_score":     integer,   // --mkini default: -1
  "gap_open_penalty":   integer,   // --mkini default: -2
  "gap_extend_penalty": integer,   // --mkini default: -1
  "max_mismatches":     integer    // --mkini default: 8, >= 0
}
```

Parameters of the `pairwise` (Smith-Waterman) aligner. The template oligo
and its reverse complement are each aligned locally against the sequence. A
hit is accepted only if it covers the whole oligo, contains no gaps, and has
at most `max_mismatches` substitutions. If both strands qualify, the hit with
fewer substitutions wins. The four scoring fields only influence which local
alignment is found.

`max_mismatches` has the same meaning as `SimpleParams.max_mismatches`: a cap
on substitutions. The two aligners differ in how they search, not in what the
cap counts. The INI keeps the two values separate, but the CLI flag
`--max-mismatches` sets both, so a file can show the same value in `pairwise`
and `simple`.

### `SimpleParams`

```jsonc
{
  "max_mismatches": integer    // --mkini default: 10, >= 0
}
```

Parameters of the bitap aligners (`simple`, `simple_simd`, `simple_cuda`). For
each sequence, the gap-free hit of the oligo with the fewest substitutions
(Hamming distance) on either strand is taken. It is accepted if it has at
most `max_mismatches` substitutions. The CLI enforces
`max_mismatches < min_oligo_length`, and `max_mismatches <= 16` for
`simple_cuda`.

### `AlignerKind`

A string enum identifying the alignment backend.

| JSON value | Description |
|---|---|
| `"pairwise"` | Smith-Waterman (rust-bio) local alignment. Omitted from the output and assumed when absent. |
| `"simple"` | Scalar bitap aligner, oligos up to 64 bp. The `--mkini` default. |
| `"simple_simd"` | AVX2-vectorized bitap. Same results as `"simple"`. |
| `"simple_cuda"` | GPU bitap. Same results as `"simple"`. Recognized on all builds for round-tripping, but only usable in builds compiled with the `cuda` feature; the CLI rejects it at startup on non-cuda builds. |

### `ThreadCount`

An externally tagged enum. Exactly one of:

| Variant | JSON representation | Description |
|---|---|---|
| Fixed | `{"Fixed": <integer>}` | A specific number of threads. The diffalign CLI always emits this variant. It resolves the configured percentage of cores (`--mkini` default 100%) to a concrete count, so the result file records the actual parallelism used. The value is machine-dependent. |
| Auto | `"Auto"` | Compatibility-only. The diffalign CLI never emits `"Auto"`; it is accepted on read so files produced by the original Tauri tool continue to parse. |

---

## `LengthResult`

One entry per processed oligo length.

```jsonc
{
  "oligo_length": integer,          // the oligo length (matches the parent key)
  "positions":    [PositionResult, ...]
}
```

| Field | Type | Description |
|---|---|---|
| `oligo_length` | integer | The oligo length for this result set. Redundant with the parent object key. |
| `positions` | array of `PositionResult` | One entry for every position `0, r, 2r, …` up to `template_length - oligo_length` (`r` = `resolution`), in ascending order. In anchored mode positions stop at `min(template_length - oligo_length, template_length - A)`, where `A` is the anchor length. This means lengths shorter than the anchor end earlier than they would without anchoring, and the array can be empty. |

---

## `PositionResult`

One entry per analyzed position within a given oligo length.

```jsonc
{
  "position":        integer,                   // 0-indexed
  "variants_needed": integer,                   // >= 0
  "analysis":        WindowAnalysisResult,
  "exclusivity":     ExclusivityResult | null
}
```

| Field | Type | Description |
|---|---|---|
| `position` | integer | Zero-indexed start position of the window on the template sequence. |
| `variants_needed` | integer | Always equal to `analysis.variants_for_threshold`; kept for compatibility. |
| `analysis` | `WindowAnalysisResult` | Variant analysis of the references at this window. |
| `exclusivity` | `ExclusivityResult` or `null` | Differential analysis of the template oligo at this window. Non-null on **every** position, including skipped ones, when differential mode is enabled; `null` otherwise. |

---

## `WindowAnalysisResult`

The variant analysis for one position/length window.

```jsonc
{
  "variants":                [Variant, ...],
  "total_sequences":         integer,    // >= 1
  "sequences_analyzed":      integer,    // >= 0
  "no_match_count":          integer,    // >= 0
  "variants_for_threshold":  integer,    // >= 0
  "coverage_at_threshold":   number,     // 0.0 - 100.0
  "skipped":                 boolean,
  "skip_reason":             string | null
}
```

| Field | Type | Description |
|---|---|---|
| `variants` | array of `Variant` | The variants for this window, in the order described under [Variant order](#variant-order). At most `var_limit` entries. |
| `total_sequences` | integer | Number of reference sequences. Always equal to the top-level `total_sequences`. |
| `sequences_analyzed` | integer | References represented by an entry of `variants`: the references that matched, minus those whose variants were dropped by `var_limit`. Equal to the sum of `variants[].count`. |
| `no_match_count` | integer | `total_sequences - sequences_analyzed`. This covers references with no accepted hit **plus** references whose variant was dropped by `var_limit`; the file does not distinguish the two. |
| `variants_for_threshold` | integer | Number of leading entries of `variants` whose percentages add up to at least `coverage_threshold`. If the whole array does not reach the threshold, this equals the length of `variants`. `0` when skipped. |
| `coverage_at_threshold` | number | Sum of the percentages of those leading variants. Below `coverage_threshold` when the threshold was not reached, for example because too many references did not match or variants were dropped by `var_limit`. `0.0` when skipped. |
| `skipped` | boolean | `true` when no reference matched this window (`sequences_analyzed == 0`). |
| `skip_reason` | string or `null` | Human-readable reason when `skipped` is `true`, otherwise `null`. Currently always `"No valid matches found in any reference sequence"`; treat it as free text. |

### Coverage calculation details

- Each variant's `percentage` is `count / total_sequences * 100`. The
  denominator is **all** references, including no-matches, so if some
  references don't match, even full coverage of the matched ones stays below
  100%.
- Percentages are accumulated **in array order**, which is not necessarily
  by descending count (see [Variant order](#variant-order)).
  `variants_for_threshold` is the first N for which the sum of the first N
  percentages is `>= coverage_threshold`, and `coverage_at_threshold` is that
  sum.
- The threshold was reached exactly when
  `coverage_at_threshold >= coverage_threshold`. The values are IEEE-754
  doubles, so this comparison is exact on the stored values, but re-adding
  the percentages yourself can differ in the last digits.

---

## `Variant`

A single variant: either a distinct matched fragment or a consensus covering
several fragments.

```jsonc
{
  "sequence":   string,
  "count":      integer,   // >= 1
  "percentage": number     // 0.0 - 100.0
}
```

| Field | Type | Description |
|---|---|---|
| `sequence` | string | Length equals the parent `oligo_length`, in template orientation (see [Conventions](#conventions)). Contains uppercase `A`, `C`, `G`, `T`, plus IUPAC ambiguity codes where a consensus method merged differing fragments. Reference sequences are not sanitized: IUPAC codes and the gap character `-` in a reference count as mismatches during the search and are copied verbatim into the fragment, so they can appear under any method. |
| `count` | integer | Number of references assigned to this variant. For `NoAmbiguities`: references whose fragment is exactly this sequence. For `FixedAmbiguities` / `Incremental`: references *first* covered by this consensus. A consensus may also match fragments already assigned to an earlier variant, but those are not counted again. Counts never overlap and add up to `sequences_analyzed`. |
| `percentage` | number | `count / total_sequences * 100`. |

Edge case: under `Incremental`, if a matched fragment itself contains
ambiguity codes copied from a reference, a fallback step can remove more
fragments than it counts. The counts can then add up to less than
`sequences_analyzed`.

### Variant order

- `NoAmbiguities`: sorted by `count` descending; ties are broken by
  `sequence` ascending.
- `FixedAmbiguities`: the order in which the greedy search selected the
  consensus sequences. This is usually, but not guaranteed to be, descending
  by `count`.
- `Incremental`: the order of selection. Counts can **increase** along the
  array, because a later step may need an ambiguity code and then cover many
  more fragments.

Keep this order when displaying the data alongside `variants_for_threshold`
and `coverage_at_threshold`, because both are computed over it. Real output
for `{"Incremental": [30, 2]}` on 10 references:

```json
"variants": [
  { "sequence": "ACGTTGCAAGCTTCGGATCC", "count": 3, "percentage": 30.0 },
  { "sequence": "ACGTTHCAAGCTTCGGATCC", "count": 6, "percentage": 60.0 },
  { "sequence": "ACGTTGCAAACTTCGGATCC", "count": 1, "percentage": 10.0 }
]
```

### IUPAC ambiguity codes

Codes that a consensus method can produce. `N` is only produced when
`exclude_n` is `false`, though it can still be copied from a reference.

| Code | Represents | Meaning |
|---|---|---|
| `R` | A or G | Purine |
| `Y` | C or T | Pyrimidine |
| `S` | G or C | Strong (3 H-bonds) |
| `W` | A or T | Weak (2 H-bonds) |
| `K` | G or T | Keto |
| `M` | A or C | Amino |
| `B` | C, G, or T | Not A |
| `D` | A, G, or T | Not C |
| `H` | A, C, or T | Not G |
| `V` | A, C, or G | Not T |
| `N` | A, C, G, or T | Any base |

---

## `ExclusivityResult`

Present when `differential_enabled` is `true`. Describes how closely the
**template oligo** at this window (`template_sequence[position ..
position + oligo_length]`) matches each exclusivity (off-target) sequence.
The variants are not tested against the exclusivity set.

Each exclusivity sequence is searched with the same aligner, strands and
acceptance rules as the references. Its mismatch count is the number of
substitutions in its best accepted hit. A sequence without an accepted hit
counts as **no-match**: it has more than `max_mismatches` substitutions, or,
under `pairwise`, its best alignment is gapped or partial. A no-match is the
desirable outcome, since the oligo is then specific against that sequence.
The cap is `simple.max_mismatches` for bitap aligners and
`pairwise.max_mismatches` for `pairwise`.

```jsonc
{
  "total_sequences":     integer,              // >= 1
  "no_match_count":      integer,              // >= 0
  "mismatch_histogram":  [MismatchBucket, ...],
  "min_mismatches":      integer | null         // null = all sequences are no-match
}
```

| Field | Type | Description |
|---|---|---|
| `total_sequences` | integer | Number of exclusivity sequences. Equal to `exclusivity_sequence_count`. |
| `no_match_count` | integer | Exclusivity sequences without an accepted hit. Equal to the `count` of the no-match bucket (`0` if there is none). |
| `mismatch_histogram` | array of `MismatchBucket` | One bucket per observed mismatch count, **sorted ascending by `mismatches`**, followed by the no-match bucket if `no_match_count > 0`. Bucket counts add up to `total_sequences`. |
| `min_mismatches` | integer or `null` | Lowest mismatch count among exclusivity sequences with an accepted hit. Lower means less specific; `0` means an off-target contains the oligo exactly. `null` means every exclusivity sequence is a no-match (fully specific — ideal). |

---

## `MismatchBucket`

A single bin in the mismatch histogram.

```jsonc
{
  "mismatches":    integer,
  "count":         integer,    // >= 1
  "example_name":  string
}
```

| Field | Type | Description |
|---|---|---|
| `mismatches` | integer | Number of substitutions, from `0` (exact match, the worst case for exclusivity) up to the mismatch cap. The sentinel value `4294967295` (u32 max) marks the no-match bucket. It is always the last bucket and fits in a JavaScript number. |
| `count` | integer | Number of exclusivity sequences in this bucket. |
| `example_name` | string | FASTA header text (after `>`) of the first exclusivity sequence, in input order, that falls in this bucket. |

---

## `Annotation`

A user annotation on the template sequence. The diffalign CLI always writes
an empty `annotations` array; this type is carried through unchanged so
files produced by the original Tauri program continue to parse and re-emit
without loss.

```jsonc
{
  "name":      string,
  "start":     integer,   // 0-indexed, inclusive
  "end":       integer,   // 0-indexed, exclusive (>= start)
  "direction": string
}
```

| Field | Type | Description |
|---|---|---|
| `name` | string | Free-form annotation label. |
| `start` | integer | Start position on the template (0-indexed, inclusive). |
| `end` | integer | End position on the template (0-indexed, exclusive). |
| `direction` | string | Free-form orientation tag (e.g. `"forward"`, `"reverse"`); diffalign does not interpret the value. |

---

## Invariants and Validation Rules

Every file written by the diffalign CLI satisfies all of the following:

1. **`template_sequence`** contains only the characters `A`, `C`, `G`, `T`
   (uppercase).
2. **`template_length`** equals the length of `template_sequence`.
3. `1 <= min_oligo_length <= max_oligo_length`.
4. The keys of `results_by_length` are exactly `min_oligo_length`,
   `min_oligo_length + s`, `min_oligo_length + 2s`, … `<= max_oligo_length`
   with `s = length_skip + 1`. Each key equals the `oligo_length` of its
   value.
5. Within a `LengthResult`, positions are exactly the multiples of
   `resolution` from `0` to `P` in ascending order, where
   `P = template_length - oligo_length`. In anchored mode,
   `P = min(template_length - oligo_length, template_length - A)` with
   `A = anchored_length` (or `min_oligo_length` when absent), and the array
   is empty when `P < 0`.
6. In each `WindowAnalysisResult`, `total_sequences` equals the top-level
   `total_sequences`, and `sequences_analyzed + no_match_count ==
   total_sequences`.
7. The sum of `variants[].count` equals `sequences_analyzed`, except in the
   `Incremental` edge case described under [`Variant`](#variant).
8. When `var_limit` is set, each `variants` array has at most `var_limit`
   entries.
9. Each variant `sequence` has length `oligo_length`, and each `percentage`
   equals `count / total_sequences * 100`.
10. `variants_needed == analysis.variants_for_threshold <= len(variants)`.
11. `skipped` is `true` exactly when `sequences_analyzed == 0`. Then
    `variants` is empty, `variants_for_threshold` is `0`,
    `coverage_at_threshold` is `0.0` and `skip_reason` is a string.
    Otherwise `skip_reason` is `null`.
12. With `NoAmbiguities`, variants are sorted by `count` descending, then
    `sequence` ascending. Other methods give no ordering guarantee.
13. `coverage_threshold` is in the range `[0.0, 100.0]`.
14. When `differential_enabled` is `false`, `exclusivity_sequence_count` is
    `null` and every `exclusivity` is `null`.
15. When `differential_enabled` is `true`, `exclusivity_sequence_count` is
    `>= 1`, and every `exclusivity` is non-null with `total_sequences` equal
    to `exclusivity_sequence_count`.
16. In each `mismatch_histogram`: buckets are sorted ascending by
    `mismatches`. Non-sentinel values are `<=` the mismatch cap. The
    `4294967295` bucket appears only if `no_match_count > 0`, and its `count`
    equals `no_match_count`. Bucket counts add up to `total_sequences`.
    `min_mismatches` is the smallest non-sentinel `mismatches`, or `null` if
    there is none.
17. `simple` is present exactly when `aligner` is `"simple"`,
    `"simple_simd"` or `"simple_cuda"`.
18. `anchored_length` is only present when `anchored` is `true`, and then
    lies within `[min_oligo_length, max_oligo_length]`.

---

## Examples

All examples below are real output from diffalign 1.0.0, run on a 30 bp
template (`ACGTTGCAAGCTTCGGATCCGTAGCTAGGC`) and small synthetic reference
sets. Only the abbreviated anchored example is shortened.

### Pairwise aligner, no omittable fields

10 references, one oligo length (20) and `resolution` 5, so the positions
are 0, 5 and 10. The aligner, `max_seeds` and `var_limit` were set
explicitly so that every omittable field is absent. This is **not** the
`--mkini` default configuration (see the next example).

- One reference lacks the first 8 bases of the template, so it is a no-match
  at positions 0 and 5 and those percentages add up to only 90%.
- One reference is stored as the reverse complement. It is counted in the
  first variant, because fragments are reported in template orientation.

```json
{
  "params": {
    "method": "NoAmbiguities",
    "pairwise": {
      "match_score": 2,
      "mismatch_score": -1,
      "gap_open_penalty": -2,
      "gap_extend_penalty": -1,
      "max_mismatches": 8
    },
    "exclude_n": true,
    "min_oligo_length": 20,
    "max_oligo_length": 20,
    "resolution": 5,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 8
    }
  },
  "template_length": 30,
  "total_sequences": 10,
  "template_sequence": "ACGTTGCAAGCTTCGGATCCGTAGCTAGGC",
  "results_by_length": {
    "20": {
      "oligo_length": 20,
      "positions": [
        {
          "position": 0,
          "variants_needed": 3,
          "analysis": {
            "variants": [
              {
                "sequence": "ACGTTGCAAGCTTCGGATCC",
                "count": 6,
                "percentage": 60.0
              },
              {
                "sequence": "ACATTGCAAGCTTCGGATCC",
                "count": 2,
                "percentage": 20.0
              },
              {
                "sequence": "ACGTTGCAAGCTTCGGACCC",
                "count": 1,
                "percentage": 10.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 9,
            "no_match_count": 1,
            "variants_for_threshold": 3,
            "coverage_at_threshold": 90.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": null
        },
        {
          "position": 5,
          "variants_needed": 2,
          "analysis": {
            "variants": [
              {
                "sequence": "GCAAGCTTCGGATCCGTAGC",
                "count": 8,
                "percentage": 80.0
              },
              {
                "sequence": "GCAAGCTTCGGACCCGTAGC",
                "count": 1,
                "percentage": 10.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 9,
            "no_match_count": 1,
            "variants_for_threshold": 2,
            "coverage_at_threshold": 90.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": null
        },
        {
          "position": 10,
          "variants_needed": 1,
          "analysis": {
            "variants": [
              {
                "sequence": "CTTCGGATCCGTAGCTAGGC",
                "count": 9,
                "percentage": 90.0
              },
              {
                "sequence": "CTTCGGACCCGTAGCTAGGC",
                "count": 1,
                "percentage": 10.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 10,
            "no_match_count": 0,
            "variants_for_threshold": 1,
            "coverage_at_threshold": 90.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": null
        }
      ]
    }
  },
  "differential_enabled": false,
  "exclusivity_sequence_count": null,
  "annotations": []
}
```

### `params` of a default configuration

The `params` object of a run that used an unmodified `--mkini` INI and no
CLI overrides. `aligner`, `simple`, `var_limit` and `max_seeds` are written
because their values differ from the omission values. `length_skip`,
`anchored` and `anchored_length` are omitted. `thread_count` depends on the
machine.

```json
"params": {
  "method": "NoAmbiguities",
  "pairwise": {
    "match_score": 2,
    "mismatch_score": -1,
    "gap_open_penalty": -2,
    "gap_extend_penalty": -1,
    "max_mismatches": 8
  },
  "aligner": "simple",
  "simple": {
    "max_mismatches": 10
  },
  "exclude_n": true,
  "min_oligo_length": 16,
  "max_oligo_length": 30,
  "resolution": 1,
  "coverage_threshold": 95.0,
  "thread_count": {
    "Fixed": 16
  },
  "var_limit": 15,
  "max_seeds": 1000
}
```

### Differential analysis, bitap aligner, `FixedAmbiguities`, `var_limit`

10 references, 5 exclusivity sequences, length 20, `resolution` 10 (positions
0 and 10). The command line was `--aligner simple --max-mismatches 4
--method fixed --fixed-ambiguities 1 --var-limit 2 --max-seeds 50`.

- `--max-mismatches` set both `pairwise.max_mismatches` and
  `simple.max_mismatches` to 4.
- The first variant is a consensus with one ambiguity code (`Y`) covering
  7 references.
- A third variant (count 1) was dropped by `var_limit = 2`. All 10
  references matched, yet `sequences_analyzed` is 9 and `no_match_count`
  is 1.
- The histograms describe the template oligo. Every value is within the cap
  of 4; two unrelated off-targets land in the no-match bucket. At position
  10, `offtarget_3` contains the oligo exactly, so `min_mismatches` is 0.

```json
{
  "params": {
    "method": {
      "FixedAmbiguities": 1
    },
    "pairwise": {
      "match_score": 2,
      "mismatch_score": -1,
      "gap_open_penalty": -2,
      "gap_extend_penalty": -1,
      "max_mismatches": 4
    },
    "aligner": "simple",
    "simple": {
      "max_mismatches": 4
    },
    "exclude_n": true,
    "min_oligo_length": 20,
    "max_oligo_length": 20,
    "resolution": 10,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 8
    },
    "var_limit": 2
  },
  "template_length": 30,
  "total_sequences": 10,
  "template_sequence": "ACGTTGCAAGCTTCGGATCCGTAGCTAGGC",
  "results_by_length": {
    "20": {
      "oligo_length": 20,
      "positions": [
        {
          "position": 0,
          "variants_needed": 2,
          "analysis": {
            "variants": [
              {
                "sequence": "ACGTTGCAAGCTYCGGATCC",
                "count": 7,
                "percentage": 70.0
              },
              {
                "sequence": "ACGTTGCAAGCTTCGAATCC",
                "count": 2,
                "percentage": 20.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 9,
            "no_match_count": 1,
            "variants_for_threshold": 2,
            "coverage_at_threshold": 90.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": {
            "total_sequences": 5,
            "no_match_count": 2,
            "mismatch_histogram": [
              {
                "mismatches": 2,
                "count": 2,
                "example_name": "offtarget_2"
              },
              {
                "mismatches": 3,
                "count": 1,
                "example_name": "offtarget_1"
              },
              {
                "mismatches": 4294967295,
                "count": 2,
                "example_name": "offtarget_4"
              }
            ],
            "min_mismatches": 2
          }
        },
        {
          "position": 10,
          "variants_needed": 2,
          "analysis": {
            "variants": [
              {
                "sequence": "CTYCGGATCCGTAGCTAGGC",
                "count": 7,
                "percentage": 70.0
              },
              {
                "sequence": "CTTCGAATCCGTAGCTAGGC",
                "count": 2,
                "percentage": 20.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 9,
            "no_match_count": 1,
            "variants_for_threshold": 2,
            "coverage_at_threshold": 90.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": {
            "total_sequences": 5,
            "no_match_count": 2,
            "mismatch_histogram": [
              {
                "mismatches": 0,
                "count": 1,
                "example_name": "offtarget_3"
              },
              {
                "mismatches": 1,
                "count": 1,
                "example_name": "offtarget_2"
              },
              {
                "mismatches": 3,
                "count": 1,
                "example_name": "offtarget_1"
              },
              {
                "mismatches": 4294967295,
                "count": 2,
                "example_name": "offtarget_4"
              }
            ],
            "min_mismatches": 0
          }
        }
      ]
    }
  },
  "differential_enabled": true,
  "exclusivity_sequence_count": 5,
  "annotations": []
}
```

### Anchored mode and `Incremental` method (abbreviated)

Anchored screening at length 22 over the range 18–25 with `length_skip` 1,
so lengths 18, 20, 22 and 24 are processed. The method is `Incremental` with
a per-step target of 95% and an ambiguity cap of 2. The `params` object is
shown in full. `results_by_length` is reduced to its position lists, which
show the anchored position range: with a 30 bp template and anchor length
22, the positions stop at 8 for every length up to 22. Without anchoring,
length 18 would run to position 12.

```jsonc
{
  "params": {
    "method": {
      "Incremental": [
        95,
        2
      ]
    },
    "pairwise": {
      "match_score": 2,
      "mismatch_score": -1,
      "gap_open_penalty": -2,
      "gap_extend_penalty": -1,
      "max_mismatches": 6
    },
    "aligner": "simple",
    "simple": {
      "max_mismatches": 6
    },
    "exclude_n": true,
    "min_oligo_length": 18,
    "max_oligo_length": 25,
    "resolution": 1,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 8
    },
    "length_skip": 1,
    "max_seeds": 100,
    "anchored": true,
    "anchored_length": 22
  },
  "template_length": 30,
  "total_sequences": 10,
  "template_sequence": "ACGTTGCAAGCTTCGGATCCGTAGCTAGGC",
  "results_by_length": {
    "18": { "oligo_length": 18, "positions": [ /* positions 0..=8 */ ] },
    "20": { "oligo_length": 20, "positions": [ /* positions 0..=8 */ ] },
    "22": { "oligo_length": 22, "positions": [ /* positions 0..=8 */ ] },
    "24": { "oligo_length": 24, "positions": [ /* positions 0..=6 */ ] }
  },
  "differential_enabled": false,
  "exclusivity_sequence_count": null,
  "annotations": []
}
```

### Skipped position

When no reference matches a window, the position is still emitted, marked as
skipped. Here 4 references cover only the first 22 bases of the template, so
the length-20 window at position 10 cannot match any of them:

```json
{
  "position": 10,
  "variants_needed": 0,
  "analysis": {
    "variants": [],
    "total_sequences": 4,
    "sequences_analyzed": 0,
    "no_match_count": 4,
    "variants_for_threshold": 0,
    "coverage_at_threshold": 0.0,
    "skipped": true,
    "skip_reason": "No valid matches found in any reference sequence"
  },
  "exclusivity": null
}
```
