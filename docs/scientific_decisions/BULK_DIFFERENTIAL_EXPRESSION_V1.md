# Bulk Treatment-versus-Matched-Control Differential Expression v1 (Step 7B)

## 1. Scientific Objective

The primary objective of Step 7B (Bulk Treatment-versus-Matched-Control Differential Expression v1) is to convert an inferentially eligible molecular-response contrast (Step 7A) into a single, standardized, gene-level treatment-response profile using official Bioconductor `DESeq2`.

For each validated treatment condition and its explicit scientifically matched controls within an approved Normalization Cohort, this stage evaluates the differential expression of all uniquely mapped canonical human genes. The resulting profile captures:
- `canonical_gene_id`: Ensembl 112 stable gene identifier.
- `approved_symbol`: HGNC approved gene symbol (when mapped).
- `base_mean`: Mean of normalized counts across all evaluated contrast samples.
- `log2_fold_change`: Unshrunk maximum likelihood estimate (MLE) of $\log_2(\text{Treatment} / \text{Control})$.
- `lfc_standard_error`: Standard error of the $\log_2\text{FC}$ estimate ($SE$).
- `wald_statistic`: Wald test statistic ($z = \log_2\text{FC} / SE$).
- `p_value`: Two-sided Wald test asymptotic p-value.
- `adjusted_p_value_bh`: Benjamini-Hochberg (BH) false discovery rate (FDR).
- `result_status`: Explicit classification (`OK`, `OUTLIER`, `ZERO_COUNTS`, `FIT_ERROR`).
- `result_notes`: Diagnostic notes explaining undefined or outlier statistics.

---

## 2. Requirement for Explicit Matched Controls

Differential expression inference requires biologically and experimentally matched baseline controls. Biological systems such as human neural organoids exhibit baseline variation across differentiation protocols, parental biological sources, culture media batches, and developmental ages.

Step 7B strictly rejects pseudo-controls, cross-experiment controls, and historical uncalibrated controls. All contrasts entering Step 7B must originate from Step 7A `MolecularResponseContrast` instances established within the same `NormalizationCohort`, sharing:
- Identical study and experiment boundaries;
- Identical bulk RNA-seq sequencing modality and library context;
- Compatible developmental age and organoid differentiation stage;
- Explicit `treatment_control_relationship_id` metadata.

---

## 3. Mathematical Directionality and Log2FC Orientation

Directionality is strictly governed:
$$\log_2\text{FoldChange} = \log_2\left(\frac{\text{Treatment}}{\text{Control}}\right)$$

- **Positive $\log_2\text{FC} > 0$** strictly indicates **Treatment $>$ Control** (upregulation in treatment).
- **Negative $\log_2\text{FC} < 0$** strictly indicates **Treatment $<$ Control** (downregulation in treatment).

To ensure this orientation without ambiguity:
1. In `DESeq2::DESeqDataSetFromMatrix`, the `condition` factor is explicitly configured with `c("control", "treatment")`, establishing `"control"` as the reference level.
2. Result extraction is executed with `DESeq2::results(dds, contrast = c("condition", "treatment", "control"), independentFiltering = FALSE)`.
3. Positive and negative control synthetic tests explicitly verify this orientation.

---

## 4. Raw Integer Counts Requirement

Step 7B operates exclusively on the preserved raw integer count matrix (`CountMatrix`).
- Diagnostic normalized counts (e.g. median-of-ratios normalized floating-point counts) must **never** be passed into DESeq2.
- TPM, FPKM, CPM, log-transformed values, and variance-stabilized values are strictly forbidden as input.
- The raw count matrix remains immutable; Step 7B constructs an isolated, derived DESeq2 input subset.

---

## 5. Reuse of Locked Size Factors (Step 6C)

Bulk Normalization v1 (Step 6C) has already estimated official DESeq2 size factors across the full approved Normalization Cohort.
- Size factors must **never** be re-estimated inside Step 7B on the contrast subset.
- In R: `sizeFactors(dds) <- supplied_locked_size_factors`.
- Every sample in the contrast must possess:
  1. Exactly one size factor;
  2. A finite floating-point value;
  3. A strictly positive value ($> 0$);
  4. Exact 1:1 sample identity correspondence.
- Any mismatch blocks execution immediately.

---

## 6. Biological Replication Definition

Valid differential expression inference requires empirical measurement of biological variance within treatment and control groups:
- **Treatment Arm:** Minimum 2 distinct biological replicates ($\ge 2$).
- **Control Arm:** Minimum 2 distinct biological replicates ($\ge 2$).
- If either arm has $< 2$ distinct biological replicates, the contrast is inferentially unqualified and blocked (`insufficient_biological_replicates`).

---

## 7. Paired vs. Unpaired Biological Source Handling (DesignResolver)

Biological source (donor, cell line, genetic background) structure is analyzed deterministically:

- **Case A: Single Biological Source / Line**
  All treatment and control replicates originate from the same biological source.
  $$\sim \text{condition}$$
  Rationale: Single biological background; condition is the sole varying factor.

- **Case B: Multiple Paired / Blocked Biological Sources**
  Multiple sources are present, and every source contributes both treatment and control samples.
  $$\sim \text{biological\_source} + \text{condition}$$
  Rationale: Biological source is included as an additive blocking covariate, provided the design matrix is full rank and residual degrees of freedom remain $> 0$.

- **Case C: Independent Unpaired Biological Sources**
  Each sample represents an independent biological replicate from a distinct donor/source without reuse (sources are unique to each sample).
  $$\sim \text{condition}$$
  Rationale: Individual biological sources serve as independent replicate units. Adding a fixed effect per source would saturate the model ($N = p$); hence sources are not added as fixed covariates.

- **Case D: Treatment/Source Confounding**
  Treatment condition is completely inseparable from biological source (e.g. all treatment replicates from Source A, all control replicates from Source B).
  **Outcome:** `BLOCKED` with `treatment_source_confounding`. DESeq2 is never run because treatment effect cannot be mathematically distinguished from background source variation.

- **Case E: Partially Crossed / Ambiguous Biological Sources**
  Sources partially overlap across arms (e.g. Source A is present in treatment and control, but Source B is present only in treatment).
  **Outcome:** `NEEDS_REVIEW` (`partially_crossed_biological_sources`). The system avoids guessing and requests scientific adjudication.

---

## 8. Technical Replicate Policy

Technical replicates (multiple sequencing libraries, flow cells, or assay runs from the same biological replicate) must **never** be passed to DESeq2 as independent biological observations (pseudoreplication).
- If multiple sample columns in the count matrix share the same `biological_replicate_id`:
  - They are not counted independently toward biological $N$.
  - They are not silently aggregated in v1 unless an upstream approved aggregation rule exists.
  - The contrast triggers `NEEDS_REVIEW` (`unresolved_technical_replicates`).

---

## 9. Estimability and Model Matrix Validation

Prior to invoking DESeq2, `DesignResolver` algebraically checks:
1. Model matrix full column rank ($\text{rank}(X) = p$);
2. Linear independence of treatment condition ($\text{rank}(X_{-\text{cond}}) < \text{rank}(X)$);
3. Residual degrees of freedom ($N - \text{rank}(X) > 0$).
Failure to satisfy any condition blocks execution before calling R.

---

## 10. Gene Harmonization Boundary and Evaluated Universe

To guarantee cross-experiment comparability and prevent false mappings:
- Differential expression results are restricted to features with `GeneMappingStatus.UNIQUELY_MAPPED` to an Ensembl 112 canonical gene identifier without unresolved collisions.
- Ambiguous mappings, unmapped identifiers, and unresolved colliding features are excluded from model-ready output.
- **No Count-Summing:** Colliding features are never merged by summing raw counts.
- **No Gene Pre-filtering:** Genes are **not** filtered prior to DESeq2 based on raw count thresholds, known DNT biology, or outcome expectations.

### Gene-Universe Conservation Audit
Genes excluded from model-ready canonical DE output because of Gene Harmonization v1 constraints are never silently dropped or lost. For every Step 7B execution, exact gene-universe conservation is strictly enforced:

$$\text{total\_input\_genes\_count} = \text{eligible\_canonical\_genes\_count} + \text{excluded\_genes\_count}$$

where:
$$\text{excluded\_genes\_count} = N_{\text{unmapped}} + N_{\text{ambiguous}} + N_{\text{invalid}} + N_{\text{collision}}$$

Every excluded feature is permanently recorded in an immutable `ExcludedGeneAudit` structure on the contrast result, preserving:
- `original_gene_id`: Source identifier from the raw count matrix.
- `source_index`: Zero-based row index in the raw count matrix.
- `mapping_status`: Harmonization status (`UNMAPPED`, `AMBIGUOUS`, `INVALID`, `COLLISION`).
- `canonical_gene_id`: Target Ensembl identifier if applicable.
- `exclusion_reason`: Deterministic justification for exclusion from canonical DE analysis.

The exclusion decision is completely independent of treatment outcome, $\log_2\text{FC}$, p-value, FDR, and DNT label.

---

## 11. DESeq2 Execution Parameters

The statistical model is fitted using official DESeq2:
- `test = "Wald"` (two-sided Wald test).
- `betaPrior = FALSE` (standard unshrunk maximum likelihood estimators).
- `minReplicatesForReplace = Inf` (outlier counts are never silently replaced with trimmed means).
- `independentFiltering = FALSE` (complete evaluable gene universe retained; no significance-dependent filtering).
- Outliers identified via Cook's distance maintain their calculated `baseMean` and `log2FoldChange` while `p_value` and `adjusted_p_value_bh` are preserved as `None` with `result_status = "OUTLIER"`.
- Genes with zero counts across all samples report `base_mean = 0.0` with `None` statistics and `result_status = "ZERO_COUNTS"`.

---

## 12. Multiple-Testing Correction

False discovery rate is controlled across all evaluable genes in the contrast using the Benjamini-Hochberg (BH) procedure.
- Step 7B outputs full profiles regardless of significance.
- No FDR cutoff (e.g. $FDR < 0.05$) is used to drop genes from the output dataset.

---

## 13. DNT Label Blindness

Step 7B is strictly blind to:
- Compound DNT status (positive, negative, reference, unclassified);
- Regulatory classifications (OECD, EPA, NTP, EFSA);
- Known DNT biomarker / target lists;
- Machine learning class assignments.

A test explicitly demonstrates that adding or altering DNT metadata produces numerically identical differential expression outputs.

---

## 14. Preservation of Developmental Age and Exposure Context

Organoid developmental age and chemical exposure parameters are critical biological covariates:
- `developmental_age_or_stage_at_exposure`
- `treatment_collection_age_or_stage`
- `control_collection_age_or_stage`
- `treatment_collection_age_normalized`
- `control_collection_age_normalized`
- `exposure_duration` & `exposure_duration_unit`
- `washout_or_recovery_duration` & `washout_or_recovery_duration_unit`
- Complete `treatment_exposures` tuple.

All exposure and stage context fields are preserved verbatim on each `DifferentialExpressionContrastResult`.

---

## 15. Future Treatment $\times$ Developmental-Age Interaction Modeling

Step 7B evaluates individual treatment-versus-matched-control contrasts at discrete stages (e.g. Day 30 vs Day 30 control, Day 60 vs Day 60 control).
It does **not** fit treatment $\times$ developmental-age interaction models in this stage.
Preserving discrete stage-specific responses allows future stages to model developmental stage shifts (e.g. $\Delta \text{Response} = \text{Response}_{\text{Day 60}} - \text{Response}_{\text{Day 30}}$) in a dedicated, revalidated framework.

---

## 16. Reproducibility and Version Locking

Execution is version-locked to the provisioned and validated R runtime:
- **R:** 4.4.3
- **Bioconductor:** 3.20
- **DESeq2:** 1.46.0

Strict version verification is enabled by default. Any discrepancy blocks production execution.
All execution provenance records input matrix hashes, coldata hashes, size factor hashes, script checksum, and runtime environment telemetry.

---

## 17. Known Limitations

1. **Unbalanced Complex Designs:** Multi-factorial experiments involving partially crossed cell lines or fractional factorial designs require expert review (`NEEDS_REVIEW`) rather than automatic resolution.
2. **Batch Effects:** Experimental batch effects not explicitly captured as balanced blocking covariates cannot be automatically corrected in v1.
3. **Small Sample Sizes:** While $N \ge 2$ per arm is inferentially admissible, power to detect small effect sizes remains constrained at $N=2$.
