# Oligoscreen Results File Format Specification

This document is the definitive reference for the JSON result file format produced and consumed by
Oligoscreen Desktop. A file that conforms to this specification can be loaded by the application via
its "Load Results" feature. The application itself writes files in this exact format via "Save Results".

## Overview

- **Format**: JSON
- **Encoding**: UTF-8
- **Serialization**: Rust `serde_json::to_string_pretty` (pretty-printed, 2-space indent)
- **File extension**: `.json` (convention; not enforced)
- **Default filename**: `screening_results.json`

## Top-Level Object: `ScreeningResults`

```jsonc
{
  "params":                      AnalysisParams,
  "template_length":             integer,       // >= 1
  "total_sequences":             integer,       // >= 1
  "template_sequence":           string,        // uppercase DNA: A, C, G, T only
  "results_by_length":           { "<length>": LengthResult, ... },
  "differential_enabled":        boolean,
  "exclusivity_sequence_count":  integer | null
}
```

| Field | Type | Description |
|---|---|---|
| `params` | `AnalysisParams` | The parameters that were used to generate these results. |
| `template_length` | integer | Length of `template_sequence` in base pairs. |
| `total_sequences` | integer | Number of reference sequences that were analyzed. |
| `template_sequence` | string | The full template DNA sequence. Uppercase characters `A`, `C`, `G`, `T` only. |
| `results_by_length` | object | Map from oligo length (as a **string** key, e.g. `"18"`) to its `LengthResult`. |
| `differential_enabled` | boolean | `true` if exclusivity sequences were provided and differential analysis was performed. Defaults to `false` when absent. |
| `exclusivity_sequence_count` | integer or `null` | Number of exclusivity sequences used. `null` when differential mode is disabled. Defaults to `null` when absent. |

> **Key serialization note**: `results_by_length` keys are stringified integers (e.g. `"18"`, `"25"`),
> not bare integer keys. This is because JSON object keys must be strings.

---

## `AnalysisParams`

```jsonc
{
  "method":             AnalysisMethod,
  "pairwise":           PairwiseParams,
  "exclude_n":          boolean,
  "min_oligo_length":   integer,   // >= 1
  "max_oligo_length":   integer,   // >= min_oligo_length
  "resolution":         integer,   // >= 1
  "coverage_threshold": number,    // 0.0 - 100.0
  "thread_count":       ThreadCount
}
```

| Field | Type | Default | Description |
|---|---|---|---|
| `method` | `AnalysisMethod` | `"NoAmbiguities"` | Which variant-finding algorithm was used. |
| `pairwise` | `PairwiseParams` | see below | Pairwise alignment scoring parameters. |
| `exclude_n` | boolean | `true` | Whether sequences containing ambiguous base `N` were excluded from analysis. |
| `min_oligo_length` | integer | `18` | Shortest oligo length tested. |
| `max_oligo_length` | integer | `25` | Longest oligo length tested. Must be >= `min_oligo_length`. |
| `resolution` | integer | `1` | Step size for position iteration. `1` = every position, `2` = every other, etc. |
| `coverage_threshold` | number | `90.0` | Target cumulative variant coverage percentage (0.0 - 100.0). |
| `thread_count` | `ThreadCount` | `"Auto"` | Thread configuration used for analysis. |

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

Scoring parameters for the pairwise sequence alignment algorithm. `max_mismatches` sets the upper
limit on mismatches for an alignment to be considered a valid match.

### `ThreadCount`

A tagged enum. Exactly one of:

| Variant | JSON representation | Description |
|---|---|---|
| Auto | `"Auto"` | Used all available CPU cores. |
| Fixed | `{"Fixed": <integer>}` | Used a specific number of threads. |

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
| `no_match_count` | integer | Number of reference sequences that did **not** match (no valid alignment). Equal to `total_sequences - sequences_analyzed`. |
| `variants_for_threshold` | integer | Number of top variants needed to cumulatively reach `coverage_threshold`. |
| `coverage_at_threshold` | number | Actual cumulative coverage percentage achieved with `variants_for_threshold` variants. Range 0.0 - 100.0. |
| `skipped` | boolean | `true` if this position was skipped during analysis. |
| `skip_reason` | string or `null` | Human-readable reason for skipping. `null` when `skipped` is `false`. |

### Coverage calculation details

- Variants are evaluated cumulatively from highest `count` to lowest.
- Each variant's `percentage` is computed as: `count / total_sequences * 100`. The denominator
  includes no-match sequences, so if many references don't align, even perfect coverage of aligned
  sequences won't reach 100%.
- `variants_for_threshold` is the first count N such that the sum of percentages of the top N
  variants >= `coverage_threshold`.

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

These appear in variant sequences only when `FixedAmbiguities` or `Incremental` methods are used:

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

Present only when `differential_enabled` is `true`. Describes how well the variants at a given
position discriminate against the exclusivity (off-target) sequences.

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
| `min_mismatches` | integer or `null` | Lowest mismatch count observed across all exclusivity sequences. `null` means all exclusivity sequences are no-match (fully specific -- ideal). |

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
| `mismatches` | integer | Number of mismatches in this bin. `0` = exact match (worst case for exclusivity -- the primer also matches off-target). The sentinel value `4294967295` (u32 max) represents "no match". |
| `count` | integer | Number of exclusivity sequences with this mismatch count. |
| `example_name` | string | Name/identifier of one representative sequence in this bin. |

---

## Invariants and Validation Rules

A valid result file must satisfy all of the following:

1. **`template_sequence`** contains only the characters `A`, `C`, `G`, `T` (uppercase).
2. **`template_length`** equals the length of `template_sequence`.
3. **`min_oligo_length` <= `max_oligo_length`**, and both are >= 1.
4. **`results_by_length`** keys (as integers) span a subset of `min_oligo_length` through
   `max_oligo_length`. Each key must equal the `oligo_length` field of its value.
5. For each position: `0 <= position <= template_length - oligo_length`.
6. Positions within a `LengthResult` are ordered ascending by `position` and spaced according to
   `resolution`.
7. `no_match_count == total_sequences - sequences_analyzed` in each `WindowAnalysisResult`.
8. `coverage_threshold` is in the range `[0.0, 100.0]`.
9. Variant `percentage` values are in `[0.0, 100.0]`.
10. Variants are sorted by `count` descending.
11. When `differential_enabled` is `false`, all `exclusivity` fields must be `null` and
    `exclusivity_sequence_count` must be `null`.
12. When `differential_enabled` is `true`, `exclusivity_sequence_count` must be a positive integer.
13. `mismatch_histogram` entries are sorted ascending by `mismatches`.

---

## Minimal Valid Example

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
    "thread_count": "Auto"
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
  "exclusivity_sequence_count": null
}
```

## Example with Differential Analysis

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
    "exclude_n": true,
    "min_oligo_length": 20,
    "max_oligo_length": 20,
    "resolution": 1,
    "coverage_threshold": 90.0,
    "thread_count": {
      "Fixed": 4
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
          "variants_needed": 1,
          "analysis": {
            "variants": [
              {
                "sequence": "ATGCATGCATGCATGCATGC",
                "count": 8,
                "percentage": 80.0
              },
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
  "exclusivity_sequence_count": 5
}
```

## Example with Skipped Position

When a position cannot be analyzed (e.g., no reference sequences align at that window):

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
