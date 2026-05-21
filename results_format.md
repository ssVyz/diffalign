# diffalign Results File Format Specification

This document is the definitive reference for the JSON result file format
produced and consumed by `diffalign`. A file that conforms to this
specification can be round-tripped through `serde_json` against the types
defined in `src/analysis/types.rs`.

## Overview

- **Format**: JSON
- **Encoding**: UTF-8
- **Serialization**: Rust `serde_json::to_string_pretty` (pretty-printed,
  2-space indent)
- **File extension**: `.json` (convention; not enforced)
- **Output path**: chosen by the user via `--output`; the binary has no
  default filename and refuses to overwrite an existing file.

## A note on omitted fields

Several `AnalysisParams` fields are annotated with
`#[serde(skip_serializing_if = ...)]` so that files written with default
settings stay byte-identical to files produced by older builds (and by the
original Tauri program). Each affected field is marked **(omitted when ...)**
below. On deserialization, every such field defaults to the value listed.

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
| `total_sequences` | integer | Number of reference sequences that were analyzed. |
| `template_sequence` | string | The full template DNA sequence. Uppercase characters `A`, `C`, `G`, `T` only. |
| `results_by_length` | object | Map from oligo length (as a **string** key, e.g. `"18"`) to its `LengthResult`. Keys are emitted in ascending numeric order (backed by a `BTreeMap`). |
| `differential_enabled` | boolean | `true` if exclusivity sequences were provided and differential analysis was performed. Defaults to `false` when absent. |
| `exclusivity_sequence_count` | integer or `null` | Number of exclusivity sequences used. `null` when differential mode is disabled. Defaults to `null` when absent. |
| `annotations` | array of `Annotation` | User annotations on the template. The diffalign CLI always writes an empty array; the field exists for compatibility with files produced by the original Tauri program. Defaults to `[]` when absent. |

> **Key serialization note**: `results_by_length` keys are stringified
> integers (e.g. `"18"`, `"25"`), not bare integer keys. This is because JSON
> object keys must be strings.

---

## `AnalysisParams`

```jsonc
{
  "method":             AnalysisMethod,
  "pairwise":           PairwiseParams,
  "aligner":            AlignerKind,         // omitted when "pairwise" (default)
  "simple":             SimpleParams,        // omitted when aligner is not a bitap variant
  "exclude_n":          boolean,
  "min_oligo_length":   integer,             // >= 1
  "max_oligo_length":   integer,             // >= min_oligo_length
  "resolution":         integer,             // >= 1
  "coverage_threshold": number,              // 0.0 - 100.0
  "thread_count":       ThreadCount,
  "length_skip":        integer,             // omitted when 0 (default)
  "var_limit":          integer,             // omitted when unset (unlimited)
  "max_seeds":          integer,             // omitted when 50 (default)
  "anchored":           boolean,             // omitted when false (default)
  "anchored_length":    integer              // omitted when unset
}
```

| Field | Type | Default | Description |
|---|---|---|---|
| `method` | `AnalysisMethod` | `"NoAmbiguities"` | Which variant-finding algorithm was used. |
| `pairwise` | `PairwiseParams` | see below | Smith-Waterman scoring parameters. Always emitted (even when a bitap aligner was selected) so the file documents the full configuration. |
| `aligner` | `AlignerKind` | `"pairwise"` | Alignment backend used to match references against the template window. Omitted from JSON output when set to the default `"pairwise"`. |
| `simple` | `SimpleParams` | absent | Parameters for the bitap aligner family. Present when `aligner` is `"simple"`, `"simple_simd"`, or `"simple_cuda"`. |
| `exclude_n` | boolean | `true` | Whether reference sequences containing ambiguous base `N` were excluded from analysis. |
| `min_oligo_length` | integer | `18` | Shortest oligo length tested. |
| `max_oligo_length` | integer | `25` | Longest oligo length tested. Must be >= `min_oligo_length`. |
| `resolution` | integer | `1` | Step size for template position iteration. `1` = every position, `2` = every other, etc. |
| `coverage_threshold` | number | `90.0` | Target cumulative variant coverage percentage (0.0 - 100.0). |
| `thread_count` | `ThreadCount` | `{"Fixed": N}` | Thread configuration used for analysis. The CLI always resolves to `Fixed`; see [`ThreadCount`](#threadcount). |
| `length_skip` | integer | `0` | Number of lengths to skip between processed lengths within `[min_oligo_length, max_oligo_length]`. `0` processes every length; `1` every other; etc. Omitted from output when `0`. |
| `var_limit` | integer | absent | Maximum number of variants recorded per position. Absent = unlimited. When set and exceeded, trailing variants are dropped and their counts fold into the no-match category for that position. |
| `max_seeds` | integer | `50` | Number of seed sequences the ambiguity-aware variant finder tries per greedy step. Only consulted by `FixedAmbiguities` and `Incremental`. Omitted from output when the default `50`. |
| `anchored` | boolean | `false` | When `true`, the per-position search ran once at `anchored_length` and reused those positions to derive matched fragments for every length in `[min_oligo_length, max_oligo_length]`. Omitted from output when `false`. |
| `anchored_length` | integer | absent | Length of the anchor search when `anchored` is `true`. Absent = falls back to `min_oligo_length`. Must lie within `[min_oligo_length, max_oligo_length]`. |

### `AnalysisMethod`

A tagged enum. Exactly one of the following forms:

| Variant | JSON representation | Description |
|---|---|---|
| No ambiguities | `"NoAmbiguities"` | Find all unique exact variants. No IUPAC ambiguity codes in output. |
| Fixed ambiguities | `{"FixedAmbiguities": <integer>}` | Allow up to N ambiguity codes per variant. Value is the max allowed count (>= 0). |
| Incremental | `{"Incremental": [<integer>, <integer or null>]}` | Each variant step targets the given percentage (0-100) of remaining uncovered sequences. Second element is optional max-ambiguities limit (`null` for unlimited). |

### `PairwiseParams`

```jsonc
{
  "match_score":        integer,   // default: 2
  "mismatch_score":     integer,   // default: -1
  "gap_open_penalty":   integer,   // default: -2
  "gap_extend_penalty": integer,   // default: -1
  "max_mismatches":     integer    // default: 8, >= 0
}
```

Smith-Waterman scoring parameters. The four scoring fields are SW-specific
and do not apply when a bitap aligner is selected. `max_mismatches` here is
the SW-score-rejection cap; it shares its name with `SimpleParams.max_mismatches`
but has a different meaning (see below).

### `SimpleParams`

```jsonc
{
  "max_mismatches": integer    // default: 8, >= 0
}
```

Parameters for the bitap (`simple` / `simple_simd` / `simple_cuda`) aligners.
`max_mismatches` is the only knob bitap needs and is interpreted as a strict
substitution count (no gaps).

### `AlignerKind`

A string enum identifying the alignment backend.

| JSON value | Description |
|---|---|
| `"pairwise"` | Smith-Waterman pairwise aligner. The default. Omitted from JSON output. |
| `"simple"` | Scalar bitap aligner. |
| `"simple_simd"` | SIMD-vectorized bitap aligner. |
| `"simple_cuda"` | GPU bitap aligner. Recognized on all builds for round-tripping; only usable in builds compiled with the `cuda` feature, and the CLI rejects it at startup on non-cuda builds. |

### `ThreadCount`

A tagged enum. Exactly one of:

| Variant | JSON representation | Description |
|---|---|---|
| Fixed | `{"Fixed": <integer>}` | A specific number of threads. The diffalign CLI always emits this variant — it resolves any user-supplied percentage to a concrete core count so the result file records the actual parallelism used. |
| Auto | `"Auto"` | Compatibility-only. The diffalign CLI never emits `"Auto"`; it is accepted on read so files produced by the original Tauri tool continue to parse. |

---

## `LengthResult`

One entry per oligo length that was analyzed.

```jsonc
{
  "oligo_length": integer,          // the oligo length (matches the parent key)
  "positions":    [PositionResult, ...]
}
```

| Field | Type | Description |
|---|---|---|
| `oligo_length` | integer | The oligo length for this result set. Redundant with the parent object key. |
| `positions` | array of `PositionResult` | Results for each analyzed template position, ordered by position ascending. |

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
| `position` | integer | Zero-indexed start position on the template sequence. Valid range: `0` to `template_length - oligo_length`. |
| `variants_needed` | integer | Minimum number of variants needed to reach `coverage_threshold`. |
| `analysis` | `WindowAnalysisResult` | Core variant analysis at this position. |
| `exclusivity` | `ExclusivityResult` or `null` | Differential analysis results. `null` when differential mode is disabled. Defaults to `null` when absent. |

---

## `WindowAnalysisResult`

The core analysis output for one position/length window.

```jsonc
{
  "variants":                [Variant, ...],
  "total_sequences":         integer,    // >= 0
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
| `variants` | array of `Variant` | All distinct sequence variants found, **sorted by `count` descending** (most frequent first). |
| `total_sequences` | integer | Total number of reference sequences in the input. |
| `sequences_analyzed` | integer | Number of reference sequences that produced a valid alignment match at this position. |
| `no_match_count` | integer | Number of reference sequences that did **not** match (no valid alignment). Equal to `total_sequences - sequences_analyzed`, plus any variants dropped by `var_limit`. |
| `variants_for_threshold` | integer | Number of top variants needed to cumulatively reach `coverage_threshold`. |
| `coverage_at_threshold` | number | Actual cumulative coverage percentage achieved with `variants_for_threshold` variants. Range 0.0 - 100.0. |
| `skipped` | boolean | `true` if this position was skipped during analysis. |
| `skip_reason` | string or `null` | Human-readable reason for skipping. `null` when `skipped` is `false`. |

### Coverage calculation details

- Variants are evaluated cumulatively from highest `count` to lowest.
- Each variant's `percentage` is computed as: `count / total_sequences * 100`.
  The denominator includes no-match sequences, so if many references don't
  align, even perfect coverage of aligned sequences won't reach 100%.
- `variants_for_threshold` is the first count N such that the sum of
  percentages of the top N variants >= `coverage_threshold`.

---

## `Variant`

A single distinct oligo sequence variant.

```jsonc
{
  "sequence":   string,
  "count":      integer,   // >= 1
  "percentage": number     // 0.0 - 100.0
}
```

| Field | Type | Description |
|---|---|---|
| `sequence` | string | The variant's nucleotide sequence. Contains uppercase `A`, `C`, `G`, `T`, and potentially IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`) when ambiguity-aware methods are used. Length equals the parent `oligo_length`. |
| `count` | integer | Number of reference sequences that contain this variant at the given position. |
| `percentage` | number | `count / total_sequences * 100`. |

### IUPAC ambiguity codes

These appear in variant sequences only when `FixedAmbiguities` or
`Incremental` methods are used:

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

Present only when `differential_enabled` is `true`. Describes how well the
variants at a given position discriminate against the exclusivity (off-target)
sequences.

```jsonc
{
  "total_sequences":     integer,              // >= 0
  "no_match_count":      integer,              // >= 0
  "mismatch_histogram":  [MismatchBucket, ...],
  "min_mismatches":      integer | null         // null = all sequences are no-match
}
```

| Field | Type | Description |
|---|---|---|
| `total_sequences` | integer | Total number of exclusivity sequences analyzed. |
| `no_match_count` | integer | Exclusivity sequences that produced no valid alignment. |
| `mismatch_histogram` | array of `MismatchBucket` | Distribution of mismatch counts, **sorted ascending by `mismatches`**. |
| `min_mismatches` | integer or `null` | Lowest mismatch count observed across all exclusivity sequences. `null` means all exclusivity sequences are no-match (fully specific — ideal). |

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
| `mismatches` | integer | Number of mismatches in this bin. `0` = exact match (worst case for exclusivity — the primer also matches off-target). The sentinel value `4294967295` (u32 max) represents "no match". |
| `count` | integer | Number of exclusivity sequences with this mismatch count. |
| `example_name` | string | Name/identifier of one representative sequence in this bin. |

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

A valid result file must satisfy all of the following:

1. **`template_sequence`** contains only the characters `A`, `C`, `G`, `T`
   (uppercase).
2. **`template_length`** equals the length of `template_sequence`.
3. **`min_oligo_length` <= `max_oligo_length`**, and both are >= 1.
4. **`results_by_length`** keys (as integers) span a subset of
   `min_oligo_length` through `max_oligo_length`, respecting `length_skip`.
   Each key must equal the `oligo_length` field of its value.
5. For each position: `0 <= position <= template_length - oligo_length`.
6. Positions within a `LengthResult` are ordered ascending by `position` and
   spaced according to `resolution`.
7. In each `WindowAnalysisResult`, `sequences_analyzed + no_match_count`
   equals the number of reference sequences considered at that window. When
   `var_limit` is in effect, dropped-variant counts fold into
   `no_match_count`.
8. `coverage_threshold` is in the range `[0.0, 100.0]`.
9. Variant `percentage` values are in `[0.0, 100.0]`.
10. Variants are sorted by `count` descending.
11. When `differential_enabled` is `false`, all `exclusivity` fields must be
    `null` and `exclusivity_sequence_count` must be `null`.
12. When `differential_enabled` is `true`, `exclusivity_sequence_count` must
    be a positive integer.
13. `mismatch_histogram` entries are sorted ascending by `mismatches`.
14. When `aligner` is `"simple"`, `"simple_simd"`, or `"simple_cuda"`,
    `simple` is present. When `aligner` is `"pairwise"` (or omitted),
    `simple` is absent.
15. When `anchored` is `true`, `anchored_length` (if present) lies within
    `[min_oligo_length, max_oligo_length]`.

---

## Minimal Valid Example

Default configuration: pairwise aligner, no ambiguities, no differential.
Note how `aligner`, `simple`, `length_skip`, `var_limit`, `max_seeds`,
`anchored`, and `anchored_length` are all absent (they are at their
default-and-omitted values).

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
    "resolution": 1,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 8
    }
  },
  "template_length": 30,
  "total_sequences": 10,
  "template_sequence": "ATGCATGCATGCATGCATGCATGCATGCAT",
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
                "sequence": "ATGCATGCATGCATGCATGC",
                "count": 8,
                "percentage": 80.0
              },
              {
                "sequence": "ATGCATGCATGCATGCATGT",
                "count": 2,
                "percentage": 20.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 10,
            "no_match_count": 0,
            "variants_for_threshold": 2,
            "coverage_at_threshold": 100.0,
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

## Example with Differential Analysis and Bitap Aligner

`FixedAmbiguities`, `simple_simd` aligner, differential mode on, and an
explicit `var_limit`. Note that `aligner` and `simple` are now present,
`max_seeds` is omitted because it is at its default of `50`, and
`length_skip` / `anchored` / `anchored_length` remain omitted.

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
      "max_mismatches": 8
    },
    "aligner": "simple_simd",
    "simple": {
      "max_mismatches": 4
    },
    "exclude_n": true,
    "min_oligo_length": 20,
    "max_oligo_length": 20,
    "resolution": 1,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 4
    },
    "var_limit": 100
  },
  "template_length": 30,
  "total_sequences": 10,
  "template_sequence": "ATGCATGCATGCATGCATGCATGCATGCAT",
  "results_by_length": {
    "20": {
      "oligo_length": 20,
      "positions": [
        {
          "position": 0,
          "variants_needed": 1,
          "analysis": {
            "variants": [
              {
                "sequence": "ATGCATGCATGCATGCATGY",
                "count": 10,
                "percentage": 100.0
              }
            ],
            "total_sequences": 10,
            "sequences_analyzed": 10,
            "no_match_count": 0,
            "variants_for_threshold": 1,
            "coverage_at_threshold": 100.0,
            "skipped": false,
            "skip_reason": null
          },
          "exclusivity": {
            "total_sequences": 5,
            "no_match_count": 2,
            "mismatch_histogram": [
              {
                "mismatches": 3,
                "count": 2,
                "example_name": "OFF_TARGET_SEQ_1"
              },
              {
                "mismatches": 5,
                "count": 1,
                "example_name": "OFF_TARGET_SEQ_3"
              },
              {
                "mismatches": 4294967295,
                "count": 2,
                "example_name": "OFF_TARGET_SEQ_4"
              }
            ],
            "min_mismatches": 3
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

## Example with Anchored Mode and Incremental Method

Anchored screening at length 22, range 18-25 with `length_skip` of 1
(processes 18, 20, 22, 24). `Incremental` method with a per-step coverage
target of 95% and an ambiguity cap of 2.

```json
{
  "params": {
    "method": {
      "Incremental": [95, 2]
    },
    "pairwise": {
      "match_score": 2,
      "mismatch_score": -1,
      "gap_open_penalty": -2,
      "gap_extend_penalty": -1,
      "max_mismatches": 8
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
      "Fixed": 16
    },
    "length_skip": 1,
    "max_seeds": 100,
    "anchored": true,
    "anchored_length": 22
  },
  "template_length": 1000,
  "total_sequences": 200,
  "template_sequence": "...",
  "results_by_length": {
    "...": "..."
  },
  "differential_enabled": false,
  "exclusivity_sequence_count": null,
  "annotations": []
}
```

## Example with Skipped Position

When a position cannot be analyzed (e.g., no reference sequences align at
that window):

```json
{
  "position": 5,
  "variants_needed": 0,
  "analysis": {
    "variants": [],
    "total_sequences": 10,
    "sequences_analyzed": 0,
    "no_match_count": 10,
    "variants_for_threshold": 0,
    "coverage_at_threshold": 0.0,
    "skipped": true,
    "skip_reason": "No valid matches found in any reference sequence"
  },
  "exclusivity": null
}
```
