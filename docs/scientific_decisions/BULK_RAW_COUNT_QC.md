# Bulk RNA-seq Raw-Count Quality Control

## Purpose

Step 6A determines whether a bulk RNA-seq raw-count dataset is structurally valid and technically characterized sufficiently to continue into downstream processing. It checks count-matrix integrity, reconciles samples with canonical metadata, calculates descriptive sample metrics, and summarizes biological-replicate similarity.

This QC layer evaluates data quality and integrity. It does not assess developmental neurotoxicity (DNT), infer DNT labels, or interpret an observed expression pattern as a toxicological effect.

## Position in the pipeline

```text
metadata contract
→ metadata validation
→ ingestion / format routing
→ bulk raw-count QC
→ gene harmonization
→ normalization
→ treatment-vs-control response
→ model-ready dataset
→ AI training / calibration / validation
```

Step 6A is a pre-normalization checkpoint. Passing QC means that the input met the active structural rules and that its technical characteristics were reported; it does not establish biological validity or suitability for every later analysis.

## Inputs

Step 6A supports bulk raw-count matrices in these local formats:

- `.csv`: comma-delimited;
- `.tsv` and `.txt`: tab-delimited.

The first column contains gene identifiers and the remaining columns contain sample counts. The asset must already have been registered and routed by Step 5 to `BULK_RAW_COUNT_PIPELINE` using compatible canonical assay and format metadata.

Step 6A does not use filenames or column-name patterns to infer treatment, control, compound, replicate, or other biological identities.

## Canonical metadata authority

The accepted metadata contract is authoritative for:

- study;
- experiment;
- sequencing assay;
- input data asset;
- sample;
- condition;
- biological replicate;
- technical replicate;
- treatment or control identity;
- compound;
- concentration;
- exposure duration;
- developmental stage.

The QC layer reconciles matrix sample identifiers with canonical sample, assay, and source-asset relationships. Biological details that are absent, unknown, not reported, or not applicable remain represented that way.

Filename-derived inference is prohibited because public archives and laboratory exports often encode incomplete, inconsistent, or locally meaningful labels in filenames. Treating those labels as biological truth could silently assign a sample to the wrong condition, inflate biological replication, or create a false treatment-control comparison.

## Structural QC

Structural checks are deterministic and do not depend on unapproved numerical thresholds.

### Errors leading to `FAIL`

| Rule or condition | Scientific and engineering rationale |
|---|---|
| Missing gene-identifier column | Counts cannot be associated with biological features. |
| No sample columns | There is no sample-level expression dataset to evaluate. |
| Nonnumeric count | Raw-count calculations require numeric observations. |
| Negative count | A raw sequencing count cannot be negative. |
| Fractional or non-integer-like raw count | The input claims raw counts, which must represent discrete count observations. |
| Duplicate sample identifier | A column could not be mapped unambiguously to one sample. |
| Missing count value | Missing is not equivalent to an observed zero. |
| Empty gene identifier | The row cannot be traced to a feature. |
| Empty or malformed row/matrix | Matrix dimensions and values cannot be interpreted safely. |
| Ambiguous sample mapping | The matrix identifier does not resolve uniquely in canonical metadata. |
| Assay/sample mismatch | A mapped matrix sample is not assigned to the declared assay. |
| Assay/source-asset mismatch | The declared source asset does not belong to the assay. |
| Missing assay | QC provenance and expected sample membership cannot be established. |

An error prevents a valid matrix result. Invalid or missing values are not replaced temporarily with zero for metric calculation.

### Review findings leading to `NEEDS_REVIEW`

| Rule or condition | Reason for review |
|---|---|
| Duplicate gene identifiers | The correct harmonization or aggregation policy is not yet approved. |
| Zero-only sample | The sample contains no detected count signal but remains preserved. |
| Unexpected matrix sample | The matrix sample is absent from canonical metadata. |
| Expected assay sample missing from the matrix | The assay and matrix contents differ and require reconciliation. |
| Unresolved technical-replicate ambiguity | A biological-replicate profile would require an unapproved aggregation or selection decision. |

Duplicate genes are never automatically summed, merged, or discarded. Technical replicates are never automatically aggregated or reduced to a representative sample.

### Informational observations

Zero-only genes generate `INFO` findings because their presence is measurable but does not, by itself, invalidate the dataset. Matrix dimensions are also reported through the dataset's `number_of_genes` and `number_of_samples` fields; the accepted implementation does not emit a separate matrix-dimension finding.

### Missing counts and observed zeroes

A zero is an observed numeric count. A missing cell contains no reported observation. Step 6A preserves this distinction: zeroes contribute to zero-count metrics, while missing values are errors and are not imputed.

## Sample-level QC metrics

For every structurally valid sample column, Step 6A calculates:

| Metric | Definition or interpretation |
|---|---|
| Library size | Sum of all gene counts for the sample. |
| Detected genes | Number of genes with a count greater than zero. |
| Zero-count genes | Number of genes with a count equal to zero. |
| Zero fraction | Zero-count genes divided by the total number of genes. |
| Minimum | Smallest gene count. |
| Median | Median gene count. |
| Maximum | Largest gene count. |
| `q25` | Linearly interpolated 25th percentile of gene counts. |
| `q50` | Linearly interpolated 50th percentile; equivalent to the reported median. |
| `q75` | Linearly interpolated 75th percentile of gene counts. |
| Dataset library ratio | Sample library size divided by the median library size across the dataset. |
| Condition library ratio | Sample library size divided by the median library size within its condition. |

If a ratio denominator is zero, the ratio is recorded as `null`/not evaluable. Step 6A never fabricates zero, infinity, or another replacement value.

The current result model does not calculate relative detected-gene ratios. Detected-gene counts are reported directly; adding relative detected-gene metrics would require a separately reviewed implementation change.

## Biological replicate QC

`condition_id` defines the biological comparison group. Within that condition, `biological_replicate_id` identifies independent biological replicates. A `technical_replicate_id` identifies repeated technical observations and does not increase biological sample size.

Files, columns, cells, aliquots, or technical replicates must not be counted as additional biological *n*. When technical replicates make a biological-replicate expression profile ambiguous, Step 6A records the comparison as not evaluable and produces a review finding. It does not aggregate or select among those technical replicates.

For conditions with at least two unambiguous biological-replicate profiles, Step 6A calculates pairwise descriptive Spearman correlation:

1. Counts are converted to ranks independently for each sample.
2. Tied counts receive their average rank.
3. Pearson correlation is calculated between the two rank vectors.

Constant and all-zero vectors have zero rank variance, so their correlation is explicitly not evaluable. It is never converted to zero. The result records the method as `spearman` and the implementation as `average_rank_pearson_v1`.

Replicate correlation is descriptive because `minimum_replicate_correlation` is currently `null`. It cannot automatically pass, warn, fail, or exclude a sample.

## Why Spearman is used at this stage

Raw counts are untransformed, strongly skewed, and influenced by sequencing depth. Spearman correlation gives an initial, comparatively robust description of whether samples preserve similar gene-count rankings without claiming that raw-count distances form an appropriate normalized sample space.

PCA and other formal sample-space or outlier assessments are intentionally deferred until an approved normalization or transformation, such as an appropriate variance-stabilizing or log-count procedure, has been selected and validated.

## Threshold policy

All scientifically contingent numerical thresholds are initially `null`, including:

- minimum library size;
- minimum detected genes;
- maximum zero-gene fraction;
- lower and upper relative-library-size ratios;
- minimum replicate correlation;
- outlier method and threshold.

`null` means that Step 6A calculates and reports the applicable metric but does not automatically classify or exclude a sample from that measurement. A threshold may be activated only after scientific justification, version control, sensitivity analysis, and validation for the intended datasets and use.

Structural invariants remain active even when numerical warning thresholds are `null`.

## Non-destructive QC principle

Step 6A never:

- deletes samples;
- deletes genes;
- merges duplicate genes;
- repairs counts;
- imputes missing values;
- aggregates technical replicates;
- selects a representative technical replicate;
- removes outliers;
- modifies canonical metadata;
- modifies raw count files.

QC produces measurements and findings only. Original identifiers, source data, and metadata remain immutable.

## Status model

Status follows exact severity precedence:

```text
any ERROR
→ FAIL

otherwise any REVIEW
→ NEEDS_REVIEW

otherwise any WARNING
→ PASS_WITH_WARNINGS

otherwise
→ PASS
```

`INFO` does not affect disposition.

- `FAIL`: at least one structural or referential error prevents valid QC interpretation.
- `NEEDS_REVIEW`: no error exists, but ambiguity or an unresolved scientific/technical decision requires review.
- `PASS_WITH_WARNINGS`: the dataset is structurally valid and has no review finding, but at least one configured warning rule was triggered.
- `PASS`: no error, review finding, or warning was produced under the active ruleset.

A status describes QC under a specific ruleset. It is not a toxicological classification.

## Provenance

The QC result retains or links:

- source asset ID;
- assay ID;
- study ID;
- experiment ID;
- condition ID;
- sample ID;
- biological and technical replicate identity at sample level;
- QC ruleset version;
- semantic configuration checksum;
- calculation implementation/version;
- Spearman method and implementation version.

The configuration checksum is computed from a canonical, key-sorted representation of the rules rather than incidental YAML formatting. Deterministic finding order, calculations, version identifiers, and checksums allow an auditor to establish which input, metadata relationships, code behavior, and scientific configuration produced a result.

## Fictional example

Consider a brain-organoid experiment with two conditions:

| Condition | Biological replicates |
|---|---|
| Vehicle control | `C1`, `C2`, `C3` |
| Compound X, 10 µM | `T1`, `T2`, `T3` |

Canonical metadata assigns each sample to its condition and preserves the compound, concentration, vehicle, exposure duration, developmental stage, assay, and replicate identities.

Suppose `T2` has lower library size than `T1` and `T3`, and its gene-count ranks have lower Spearman correlation with those replicates. Step 6A reports `T2`'s library metrics and pairwise correlations. Because the library-size and correlation thresholds are `null`, it does not automatically warn, exclude, or delete `T2`.

The observation is not evidence that Compound X caused developmental neurotoxicity. It could reflect technical variation, biological variation, or another factor and must be evaluated under later approved scientific procedures.

## Explicit exclusions

Step 6A does not perform:

- normalization or transformation;
- PCA;
- differential expression;
- batch correction;
- gene harmonization or duplicate-gene aggregation;
- technical-replicate aggregation;
- treatment-control response construction;
- pathway analysis;
- feature engineering;
- AI or model training;
- DNT interpretation or label assignment.

## Known limitations

- Matrices are loaded into memory.
- CSV input is comma-delimited.
- TSV and TXT inputs are tab-delimited.
- Raw-rank Spearman correlation is descriptive.
- No technical-replicate aggregation policy is approved.
- No numerical warning or outlier threshold is active.
- Duplicate-gene resolution is not implemented.
- Compressed count matrices are outside scope.
- Remote matrix retrieval is outside scope.

## Scientific decisions still pending

The following require future scientific approval:

- numerical library-size thresholds;
- detected-gene thresholds;
- zero-fraction thresholds;
- replicate-correlation threshold;
- outlier strategy;
- technical-replicate policy;
- duplicate-gene harmonization policy;
- post-normalization PCA and outlier policy.
