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

## 8. Biological Source Identifier Governance

Biological source identity can be documented across multiple standard metadata fields:
- `biological_source_id`
- `cell_line_id`
- `donor_id`
- `source_id`

`DesignResolver` and `DifferentialExpressionPipeline` adhere to strict governance rules:
1. **Deterministic Alignment:** Both modules share unified extraction logic via `resolve_sample_biological_source`.
2. **Conflict Detection:** If a sample contains multiple candidate keys with conflicting values (e.g. `donor_id != biological_source_id`), execution is blocked (`conflicting_biological_source_identifiers`).
3. **No Fabricated Fallbacks:** Fabricated identifiers (such as `"source_1"`) are strictly forbidden. If multiple biological sources exist across the cohort/contrast, any sample lacking a source identifier immediately triggers `missing_biological_source_identifier` (`Severity.ERROR`, `ContrastStatus.BLOCKED`).
4. **Immutable Assignment Binding:** The exact mapping `sample_source_assignments: dict[str, str]` evaluated by `DesignResolver` is passed directly into `coldata` for DESeq2 model fitting.

---

## 9. Experimental Batch Information & Confounding Gating

When experimental batch metadata is tracked (`experimental_batch_id`, `batch_id`, `sequencing_batch_id`, `processing_batch_id`, `library_prep_batch_id`):
1. **Conflict Detection:** Conflicting batch values on a single sample block execution (`conflicting_batch_identifiers`).
2. **Complete Tracking Requirement:** If batch is specified for some samples, all samples in the contrast must possess batch metadata; otherwise execution blocks (`missing_batch_identifier`).
3. **Treatment–Batch Confounding Check:** If treatment and control arms are completely separated across disjoint batches (e.g. all treatment in batch 1, all control in batch 2), treatment effect is mathematically indistinguishable from technical batch effects. Execution is strictly blocked (`treatment_batch_confounding`, `Severity.ERROR`, `ContrastStatus.BLOCKED`).
4. **Balanced Batch Blocking:** When batches are unconfounded and balanced across arms, `DesignResolver` includes batch as an additive blocking covariate:
   $$\sim \text{batch} + \text{condition} \quad \text{or} \quad \sim \text{biological\_source} + \text{batch} + \text{condition}$$
   provided the model matrix is full column rank ($\text{rank}(X) = p$), treatment effect is estimable, and residual degrees of freedom remain $> 0$.
5. **Rank Deficiency Safeguard:** If adding batch results in rank deficiency or eliminates residual degrees of freedom, execution is blocked (`rank_deficient_design` or `insufficient_residual_degrees_of_freedom`).

---

## 10. Technical Replicate Policy

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

## 16. Preservation of Developmental Age, Exposure Context, and Unit Canonicalization

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

### Exposure Unit Canonicalization
Units for concentration (`uM`, `µM`, `micromolar`, `mM`, `nM`, `M`, `ug/mL`, `mg/mL`, `ng/mL`) and exposure duration (`h`, `hours`, `d`, `days`, `w`, `weeks`, `min`, `minutes`) are canonically normalized to prevent spurious mismatches between synonymous notations without altering scientific meaning.

### Incompatible Exposure Conditions Safeguard
If samples within a treatment arm exhibit conflicting doses, exposure durations, or test agents without an approved sample comparison exception, the contrast is blocked (`incompatible_exposure_conditions`). Contrasts strictly represent single, coherent exposure conditions.

---

## 17. Future Treatment $\times$ Developmental-Age Interaction Modeling

Step 7B evaluates individual treatment-versus-matched-control contrasts at discrete stages (e.g. Day 30 vs Day 30 control, Day 60 vs Day 60 control).
It does **not** fit treatment $\times$ developmental-age interaction models in this stage.
Preserving discrete stage-specific responses allows future stages to model developmental stage shifts (e.g. $\Delta \text{Response} = \text{Response}_{\text{Day 60}} - \text{Response}_{\text{Day 30}}$) in a dedicated, revalidated framework.

---

## 18. Reproducibility and Version Locking

Execution is version-locked to the provisioned and validated R runtime:
- **R:** 4.4.3
- **Bioconductor:** 3.20
- **DESeq2:** 1.46.0

Strict version verification is enabled by default. Any discrepancy blocks production execution.
All execution provenance records input matrix hashes, coldata hashes, size factor hashes, script checksum, and runtime environment telemetry.

---

## 19. Real-Data Benchmark and Independent Verification

Stage 7B was validated against the real human cortical organoid dataset GSE232218 (Valproic Acid 200 µM vs Vehicle control, $N=4$ treatment vs $N=4$ control):
- **Raw Matrix Dimensions:** 60,628 input rows $\times$ 15 samples.
- **Canonical Gene Universe:** 60,207 evaluable Ensembl 112 canonical genes.
- **Excluded Features:** 421 non-canonical/unmapped features (100% audited).
- **Gene Universe Conservation:** Conserved exactly ($60,207 + 421 = 60,628$).
- **Significant Differentially Expressed Genes:** 111 upregulated, 16 downregulated (FDR $< 0.05$).
- **Independent Direct R / DESeq2 Execution Agreement:**
  - Max $\log_2\text{FC}$ difference: $4.9998 \times 10^{-9}$
  - Max Wald statistic difference: $4.9999 \times 10^{-9}$
  - Max p-value difference: $4.9995 \times 10^{-9}$
  - Max adjusted p-value difference: $4.9948 \times 10^{-9}$
  - Direction agreement: $32,227 / 32,227$ ($100.0\%$)
- **Numerical Reproducibility:** Exact agreement against verified baseline registry across all 60,207 genes ($0.000000 \times 10^0$ max diff).
- **Total Repository Test Suite:** 344 / 344 automated tests passing.

---

## 20. Known Limitations

1. **Unbalanced Complex Designs:** Multi-factorial experiments involving partially crossed cell lines or fractional factorial designs require expert review (`NEEDS_REVIEW`) rather than automatic resolution.
2. **Unbalanced Batch Effects:** Complex fractional batch designs without balance across arms require manual review.
3. **Small Sample Sizes:** While $N \ge 2$ per arm is inferentially admissible, power to detect small effect sizes remains constrained at $N=2$.
4. **Real-Data Multi-Batch Validation:** While batch-blocking formulas (`~ batch + condition` and `~ biological_source + batch + condition`) and confounding safeguards are mathematically and programmatically verified via automated tests, empirical validation on real multi-batch human organoid datasets remains to be conducted when an authenticated multi-batch dataset is ingested.

---

## 21. Researcher User Guide: Installation, Input Formats, and Verification

### 21.1 Environment Setup and Dependencies
Stage 7B requires:
- **Python:** 3.11+ (verified on 3.14.6) with packages: `pandas`, `pytest`, `pyyaml`.
- **R:** 4.4.3 with Bioconductor 3.20 and `DESeq2` 1.46.0.
- Verify R/Bioconductor setup:
  ```powershell
  Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))"
  # Output must report: 1.46.0
  ```

### 21.2 Raw Integer Counts vs. Processed Expression Routes
- **Raw Count Route (Stage 7B):** Stage 7B operates exclusively on unnormalized non-negative integer count matrices (HTSeq, featureCounts, STAR gene counts, RSEM expected counts rounded to integer). It estimates dispersion and fits negative-binomial generalized linear models.
- **Processed Expression Route (TPM, FPKM, CPM):** Processed expression matrices containing normalized, fractional, or log-transformed values cannot be modeled using DESeq2's negative-binomial distribution. They are routed upstream to the processed-expression pipeline (e.g. Stage 8A limma-trend / limma-voom) and are strictly rejected by Stage 7B.

### 21.3 Providing Sample Metadata
Sample metadata must be supplied as a CSV or TSV file containing at minimum:
- `sample_id`: Exact string matching the count matrix sample column header.
- `condition_type` (or `condition` / `treatment_control_status`): `"treatment"`, `"control"`, or `"exclude"`.
- `biological_replicate_id`: Distinct identifier per true biological replicate (e.g. `rep_trt_1`, `rep_ctrl_1`).
- `biological_source_id` (optional): Donor / cell line identifier (e.g. `WA09`, `donor_A`).
- `batch_id` (optional): Experimental / sequencing batch identifier (e.g. `batch_1`, `batch_2`).

Example template:
```csv
sample_id,condition_type,biological_replicate_id,biological_source_id,batch_id
Sample_1,control,rep_ctrl_1,line_A,batch_1
Sample_2,control,rep_ctrl_2,line_A,batch_2
Sample_3,treatment,rep_trt_1,line_A,batch_1
Sample_4,treatment,rep_trt_2,line_A,batch_2
```

### 21.4 Running Independent Direct R / DESeq2 Verification
The reusable verification engine can be executed via CLI:
```powershell
# Option A: Run using a benchmark configuration file
python scripts/independent_deseq2_verification.py --config config/benchmarks/gse232218_vpa_200um.json

# Option B: Run on custom RNA-seq count matrix and sample table
python scripts/independent_deseq2_verification.py `
  --counts path/to/raw_counts.tsv `
  --samples path/to/samples.csv `
  --agent-name "Compound_X" `
  --concentration 10.0 `
  --concentration-unit "uM" `
  --developmental-age "Day 30" `
  --exposure-duration 48.0 `
  --output-dir path/to/verification_artifacts
```

### 21.5 Gating and Status Taxonomy
Contrasts are evaluated against strict scientific invariants:

| Condition Encountered | Emitted Rule ID | Status | Rationale |
| :--- | :--- | :--- | :--- |
| Insufficient replicates ($N < 2$) | `insufficient_biological_replicates` | `BLOCKED` | Empirical variance unestimable |
| Developmental age mismatch | `developmental_age_mismatch` | `BLOCKED` | Cross-developmental comparisons invalid |
| Confounded treatment & source | `treatment_source_confounding` | `BLOCKED` | Treatment inseparable from genetics |
| Confounded treatment & batch | `treatment_batch_confounding` | `BLOCKED` | Treatment inseparable from batch |
| Conflicting source IDs | `conflicting_biological_source_identifiers` | `BLOCKED` | Inconsistent sample metadata |
| Conflicting batch IDs | `conflicting_batch_identifiers` | `BLOCKED` | Inconsistent sample metadata |
| Incompatible exposures in arm | `incompatible_exposure_conditions` | `BLOCKED` | Mixed doses/durations without exception |
| Rank deficient design matrix | `rank_deficient_design` | `BLOCKED` | Linear collinearity among covariates |
| Zero residual degrees of freedom | `insufficient_residual_degrees_of_freedom` | `BLOCKED` | Model saturated |
| Technical replicates present | `unresolved_technical_replicates` | `NEEDS_REVIEW` | Pseudoreplication avoidance |
| Partially crossed sources | `partially_crossed_biological_sources` | `NEEDS_REVIEW` | Unbalanced multi-source structure |
| Ambiguous developmental age | `unresolved_developmental_age_equivalence` | `NEEDS_REVIEW` | Missing or non-canonical age token |
