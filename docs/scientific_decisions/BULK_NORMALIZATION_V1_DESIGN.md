# Bulk RNA-seq Normalization v1: Scientific and Architectural Design

## 1. Purpose and Scope

Bulk RNA-seq Normalization v1 provides a rigorous, reproducible, provenance-preserving normalization layer for human bulk RNA-seq data that has completed bulk raw-count quality control (Step 6A) and Gene Harmonization v1 (Step 6B). Its sole objective is to estimate sequencing-depth and library-composition scaling factors (size factors, $s_j$) for samples within scientifically coherent cohorts using the official Bioconductor `DESeq2` median-of-ratios methodology, calculate secondary normalized count matrices purely for technical diagnostics and inspection, and evaluate normalization-specific technical quality control.

This stage strictly enforces that raw integer counts remain the immutable, primary modeling substrate for all downstream negative-binomial generalized linear modeling (Step 7 Treatment-vs-Control Response Construction and Differential Expression).

### Strict Pipeline Boundaries

Bulk RNA-seq Normalization v1 is strictly bounded. It does **not** perform:
- Differential expression analysis (no Wald tests, likelihood ratio tests, or dispersion shrinkage);
- Log2 fold change ($\log_2\text{FC}$) or shrinkage computation;
- P-value, FDR, or Benjamini-Hochberg adjustment calculations;
- Variance-stabilizing transformation (VST) or regularized log (rlog) transformation;
- Batch correction (no ComBat, sva, or RUVSeq);
- TPM / FPKM / CPM conversion;
- Arbitrary low-count filtering;
- Pathway, gene ontology, or functional enrichment analysis;
- DNT marker gene discovery, biomarker scoring, or phenotypic profiling;
- Feature selection, filtering, or dimensionality reduction (no PCA, UMAP, t-SNE);
- Dataset splitting;
- AI model training, calibration, validation, or inference.

---

## 2. Position in the Pipeline Architecture

```text
metadata contract
→ metadata validation
→ ingestion / format routing
→ bulk raw-count QC (Step 6A)
→ gene harmonization (Step 6B)
→ bulk RNA-seq normalization (Step 6C: Normalization v1)
  ├── 1. Normalization Cohort dynamic formation
  ├── 2. Canonical collision gating & mapping eligibility filtering
  ├── 3. Official Bioconductor DESeq2 execution (estimateSizeFactors)
  ├── 4. Primary output: sample size factors + estimation provenance
  ├── 5. Secondary output: diagnostic normalized count matrix (inspection only)
  └── 6. Treatment-control matching & ResponseEligibility tagging
→ treatment-vs-control response construction (Step 7: raw counts + size factors + design)
→ standardized feature engineering
→ model-ready dataset
→ AI training / calibration / validation
```

---

## 3. Dynamic Normalization Cohort Architecture

### 3.1 Replacement of the Rigid (Study, Experiment, Assay) Unit

In heterogeneous public transcriptomic datasets, sequencing assays (`assay_id`) frequently partition samples along operational lines (such as separate sequencing lanes, flow cells, or library batches). However, biologically matched treatment and control samples from a single planned experiment are often split across adjacent assays. Conversely, samples from completely separate studies or biological modalities might erroneously share assay identifiers.

Normalization v1 abandons the rigid `(study_id, experiment_id, assay_id)` tuple and establishes the **Normalization Cohort** ($C$). A Normalization Cohort represents a metadata-derived set of samples that are scientifically and technically appropriate for joint size-factor estimation.

### 3.2 Mandatory Cohort Invariants

The cohort constructor must strictly satisfy the following rules:

1. **Study and Experiment Boundaries:**
   - Every sample in a cohort must share identical `study_id` and `experiment_id`.
   - **Cross-study pooling is strictly prohibited.** Under no circumstances may samples from separate studies be pooled into a shared normalization cohort.
   - Experiments within the same study are not pooled unless an explicit metadata assertion establishes they belong to a single unified experimental design.
2. **Modality Homogeneity:**
   - Every sample in the cohort must have `sequencing_assay.rna_seq_modality == "bulk_rna_seq"`.
   - Single-cell (`scRNA_seq`), single-nucleus (`snRNA_seq`), or spatial transcriptomic assays must never be co-normalized with bulk RNA-seq.
3. **Sequencing and Library Context Compatibility:**
   - Samples must share compatible technical contexts based on metadata fields:
     - `library_strategy` (e.g., polyA selection vs. ribo-depletion; distinct strategies cannot share composition size factors);
     - `sequencing_platform` (e.g., Illumina NovaSeq vs. Oxford Nanopore);
     - `strandedness` (e.g., unstranded, forward, reverse);
     - `read_layout` (single-end vs. paired-end).
   - If technical context metadata is missing or conflicting across samples, the cohort constructor must **not guess**. It must flag `cohort_technical_context_ambiguous` and halt at `NEEDS_REVIEW`.
4. **Preservation of Matched Treatment and Controls:**
   - For every planned comparison defined in `treatment_control_relationships`, all treatment conditions and their corresponding matched control conditions must remain in the **same** Normalization Cohort.
   - An `assay_id` is relevant metadata but **must not automatically separate treatment and matched-control samples that belong to the same scientifically coherent normalization cohort**.
5. **Biological Context Consistency:**
   - Samples within a cohort must originate from compatible biological systems (`biological_source_id`, `organoid_context_id`).
6. **Ambiguity Gating:**
   - If available metadata is insufficient to confirm that samples can be safely co-normalized, the pipeline must emit a finding (`cohort_assignment_ambiguous`, `Severity.REVIEW`) and halt at `Status.NEEDS_REVIEW`.

### 3.3 Cohort Definition Data Model

```python
@dataclass(frozen=True, slots=True)
class NormalizationCohort:
    """Scientifically validated group of samples for joint size-factor estimation."""
    cohort_id: str
    study_id: str
    experiment_id: str
    sample_ids: tuple[str, ...]
    assay_ids: tuple[str, ...]
    rna_seq_modality: str
    library_strategy: str | None
    sequencing_platform: str | None
    strandedness: str | None
    read_layout: str | None
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]
    unassigned_sample_ids: tuple[str, ...]
    treatment_control_relationship_ids: tuple[str, ...]
    metadata_hash: str
```

---

## 4. Upstream Gating, Feature Eligibility, and Canonical Collision Gate

### 4.1 Upstream QC and Harmonization Gating

Normalization v1 requires that input data has passed upstream stages:
1. **Raw-Count QC (Step 6A):** Must be `QCStatus.PASS` or `QCStatus.PASS_WITH_WARNINGS`. An upstream `QCStatus.NEEDS_REVIEW` requires a validated `QCReviewDisposition`; `QCStatus.FAIL` blocks execution completely.
2. **Gene Harmonization (Step 6B):** Must be `Status.PASS` or `Status.PASS_WITH_WARNINGS`.

### 4.2 Strict Canonical Collision Gate

In Gene Harmonization v1, a canonical collision occurs when multiple distinct input rows map to the exact same canonical Ensembl Gene ID.
- **Default Policy: Halt at `NEEDS_REVIEW`.** There is no default automatic collision resolution.
- If an input dataset contains any canonical collisions (`len(harmonized_dataset.collisions) > 0`), normalization stops immediately with finding `canonical_collisions_detected` (`Severity.REVIEW`), dataset status `Status.NEEDS_REVIEW`.
- **Preservation of Raw Data on Collision Block:**
  - When blocked by collisions, the pipeline preserves the immutable raw/harmonized input matrix (`raw_count_matrix`);
  - Sets normalized output to unavailable/`None` (`diagnostic_normalized_matrix = None`);
  - Preserves all collision records;
  - Returns the appropriate structured review finding;
  - Does **not** represent the source raw matrix as lost.
- **Prohibition of Automatic Aggregation:** The pipeline must **never** automatically sum, average, select by mean, select by variance, choose the first row, or delete a colliding row.
- **Authorized Disposition:** Processing a dataset with canonical collisions is only permissible if an explicit, auditable, and scientifically justified **`CollisionDisposition`** is provided:
  ```python
  @dataclass(frozen=True, slots=True)
  class CollisionDisposition:
      decision: str  # e.g., "DROP_COLLISIONS_FROM_ESTIMATION", "MANUAL_REPAIR"
      rationale: str
      reviewer: str
      review_date: str
      resolved_collision_ids: tuple[str, ...]
      provenance: dict[str, Any] | None = None
  ```

### 4.3 Feature Eligibility for Size-Factor Estimation

Size factors must reflect genuine global sequencing depth and library composition. Features that are ambiguous, unmapped, or colliding introduce artificial noise into the geometric mean calculation.

- **Eligible Feature Set ($E$):**
  A row $i$ is eligible for size-factor estimation if and only if:
  1. `mapping_status == GeneMappingStatus.UNIQUELY_MAPPED`;
  2. `canonical_gene_id` is non-null and references a valid Ensembl 112 canonical locus;
  3. The row is **not** involved in any unresolved canonical collision.
- **Preservation of All Source Rows:** All other source rows and their mapping states (`AMBIGUOUS`, `UNMAPPED`, `INVALID`, colliding) remain completely preserved.
- **Strict Prohibition of Biased Feature Subsetting:**
  Do **not** use:
  - DNT labels;
  - DNT-positive / DNT-negative status;
  - Expected biomarkers;
  - The proposal marker panel;
  - Treatment outcomes;
  - Future AI features
  to choose normalization genes. Size-factor estimation must rely exclusively on objective mapping status.
- **Two Matrix Definitions:**
  1. **Estimation Submatrix ($K_{\text{eligible}}$):** Subset of rows $i \in E$, dimensions $|E| \times M$, passed directly to DESeq2 for size-factor estimation.
  2. **Full Diagnostic Matrix ($\tilde{K}$):** All original rows $i \in \{1, \dots, N\}$, scaled by the cohort size factors: $\tilde{K}_{ij} = K_{ij} / s_j$. Preserved for technical auditing and inspection only.

---

## 5. Official Bioconductor DESeq2 Execution Engine

### 5.1 Official DESeq2 Implementation Only

Production size-factor estimation must be performed by the official R/Bioconductor `DESeq2` package (`DESeq2::estimateSizeFactors(..., type="ratio")`).
- Python orchestrates validation, cohort construction, execution, serialization, provenance, and result checking.
- In-house Python reimplementations of the median-of-ratios algorithm must **not** be used as the production result.
- A controlled headless R execution interface (`src/normalization/r_bridge.py` and `scripts/r/estimate_size_factors.R`) is implemented.

### 5.2 Version-Locked R/Bioconductor Environment

The Normalization v1 reference environment is version-locked rather than specified by loose ranges:
- **R Version:** `4.4.3`
- **Bioconductor Release:** `3.20`
- **DESeq2 Version:** `1.46.0`

### 5.3 Deterministic Dependency Management and Runtime Capture

1. **Environment Definition:**
   - The environment is deterministically pinned in `scripts/r/renv.lock` and `environment-r.yaml` (specifying exact CRAN and Bioconductor repository packages and hashes).
2. **Runtime Metadata Capture:**
   At runtime, the execution bridge records the actual:
   - R version;
   - Bioconductor version;
   - DESeq2 version;
   - Relevant dependency/environment identity (e.g. `renv` lockfile hash or conda prefix);
   - Operating system;
   - Execution timestamp.
3. **Strict Version Mismatch Gating by Default:**
   - Runtime version checking is strict by default (`strict_version_check = True`).
   - A mismatch in any of the three locked versions (R `4.4.3`, Bioconductor `3.20`, DESeq2 `1.46.0`) generates a `Severity.ERROR` finding and blocks production Normalization v1 execution (`Status.FAIL`).
   - The system will **not** silently continue under an unrecorded, mismatched, or alternative software version.
   - Future software-version upgrades require an explicitly versioned and revalidated normalization contract rather than being treated as equivalent to v1.
4. **Programmatic Bioconductor Detection without Inferred Fallback:**
   - The R environment probe must programmatically verify the Bioconductor version from the runtime (via `BiocManager::version()` or `BiocVersion`).
   - The probe must **never** report Bioconductor 3.20 merely because DESeq2 is installed.
   - If the Bioconductor version cannot be reliably determined, it is reported as `unavailable`, an explicit finding (`bioc_version_undetermined`, `Severity.ERROR`) is emitted, and strict Normalization v1 environment validation fails.
5. **Environment Availability Handling:**
   - If R or DESeq2 is unavailable in the execution environment, the system does not substitute a Python approximation and does not claim scientific validation. It reports the missing runtime requirement explicitly via `r_environment_missing` (`Severity.ERROR`). Unit tests that do not require R may still execute.

### 5.4 Execution Contract and Serialization

1. **Input Serialization:**
   - Matrix written to temporary scratch directory as a clean tab-delimited file.
   - Rows: eligible canonical Ensembl Gene IDs. Columns: sample IDs in exact cohort order.
   - Counts: discrete non-negative integers.
   - SHA-256 computed on the serialized input matrix prior to R invocation.
2. **Headless R Runner (`scripts/r/estimate_size_factors.R`):**
   - Invoked via `subprocess.run(["Rscript", "--vanilla", script_path, input_tsv, method, output_json], check=False, capture_output=True)`.
   - Executes:
     ```R
     library(DESeq2)
     counts <- as.matrix(read.delim(input_tsv, row.names=1, check.names=FALSE))
     sf <- DESeq2::estimateSizeFactorsForMatrix(counts, type=method)
     ```
   - Emits structured JSON containing calculated size factors, R runtime versions, OS information, and diagnostic messages.

---

## 6. Mathematical Methodology: Standard Ratio vs. Poscounts

### 6.1 Default Method: Standard Median-of-Ratios (`type = "ratio"`)

The default normalization method is `DESeq2::estimateSizeFactors(..., type="ratio")`.
1. For eligible genes $i \in E$ across $m$ samples in cohort $C$:
   $$g_i = \left( \prod_{j=1}^m K_{ij} \right)^{1/m} = \exp\left( \frac{1}{m} \sum_{j=1}^m \ln K_{ij} \right)$$
2. Filter for genes with strictly non-zero geometric mean: $E^+ = \{i \in E \mid g_i > 0\}$.
3. Compute sample size factors:
   $$s_j = \operatorname{median}_{i \in E^+} \left( \frac{K_{ij}}{g_i} \right)$$

### 6.2 Standard Ratio Failure and Non-Automatic Fallback

If every gene has at least one sample with zero counts ($K_{ij} = 0$), $g_i = 0$ for all genes and standard DESeq2 fails with:
`"every gene contains at least one zero, cannot compute log geometric means"`.
- **Governing Behavior:**
  - Do **not** silently substitute another method.
  - Preserve the raw input (`raw_count_matrix`).
  - Produce no normalized matrix (`diagnostic_normalized_matrix = None`).
  - Record the DESeq2 failure reason in a structured finding: `standard_size_factors_failed` (`Severity.REVIEW`).
  - Place the cohort in an explicit review state (`Status.NEEDS_REVIEW`).

### 6.3 Explicit Authorization for Poscounts (`type = "poscounts"`)

`type = "poscounts"` must **never** be invoked automatically.
$$g_i^{\text{pos}} = \left( \prod_{j: K_{ij} > 0} K_{ij} \right)^{1 / m_i^+}, \quad m_i^+ = \sum_{j=1}^m \mathbb{I}(K_{ij} > 0)$$
It requires explicit, auditable scientific authorization:
- Explicit configuration `size_factor_method: poscounts` in `NormalizationConfig`; OR
- An auditable **`NormalizationReviewDisposition`** recording scientific justification, reviewer identity, and review date.

---

## 7. Primary Modeling Contract: Raw Counts vs. Size Factors

### 7.1 Preserved Output Entities

Normalization v1 preserves separately:
1. **Immutable raw integer counts:** Original harmonized counts ($K_{ij}$);
2. **DESeq2 sample size factors:** $s_j$ and their exact estimation provenance;
3. **Normalized floating-point counts:** $\tilde{K}_{ij} = K_{ij} / s_j$ for QC, inspection, and exploratory downstream use;
4. **Exact normalization-cohort membership:** Linkage to `NormalizationCohort`;
5. **Normalization method:** Record of method (`ratio` or `poscounts`);
6. **Complete execution provenance:** R/Bioconductor/DESeq2 versions, environment hash, execution timestamp, OS.

### 7.2 Downstream Differential Expression Contract

Future DESeq2 differential-expression analysis (Step 7) will consume the **appropriate raw integer counts + approved size factors/normalization contract + experimental design**.
- DESeq2 uses generalized linear modeling:
  $$K_{ij} \sim \operatorname{NB}(\mu_{ij}, \alpha_i), \quad \mu_{ij} = s_j q_{ij}, \quad \log_2(q_{ij}) = x_{j\cdot}\beta_i$$
- **Strict Prohibition:** Do **not** design later differential-expression analysis to feed normalized floating-point counts back into DESeq2 as if they were raw counts. Doing so violates discrete count distribution assumptions and corrupts dispersion estimation.

---

## 8. Missing Matched Control Semantics

### 8.1 Independence of Normalization from Control Presence

Lack of a matched control is **not** by itself a mathematical normalization failure.
- Size factors can and should still be computed across available technically valid samples in the cohort.

### 8.2 Downstream Response Eligibility Contract

Samples are tagged with an explicit **`ResponseEligibility`** status:

```python
class ResponseEligibility(StrEnum):
    ELIGIBLE_WITH_MATCHED_CONTROL = "ELIGIBLE_WITH_MATCHED_CONTROL"
    INELIGIBLE_MISSING_MATCHED_CONTROL = "INELIGIBLE_MISSING_MATCHED_CONTROL"
    CONTROL_BASELINE = "CONTROL_BASELINE"
    UNASSIGNED_OR_AMBIGUOUS = "UNASSIGNED_OR_AMBIGUOUS"
```

- `ResponseEligibility.INELIGIBLE_MISSING_MATCHED_CONTROL`:
  This flag means the sample/condition cannot contribute to the planned treatment-versus-matched-control response construction (Step 7) until an appropriate control relationship exists.
- The pipeline does not couple Normalization v1 to a generic differential expression concept; it specifically bounds treatment-vs-matched-control response construction.

---

## 9. Zero Low-Count Gene Filtering Policy

Normalization v1 performs **zero low-count gene filtering**.
- All genes present in the harmonized input must be preserved in the output.
- Any filtering based on expression level (e.g. minimum count thresholds across replicates) belongs to downstream differential expression or response construction, where experimental design and group structure can be considered properly.

---

## 10. Scientific QC Metrics & Configurable Thresholds

### 10.1 Elimination of Invented Thresholds

Calculate useful diagnostics without attaching `PASS/REVIEW/FAIL` behavior to arbitrary cutoffs (such as `< 100` genes or `> 20x` size-factor disparity):
- Number of genes used for size-factor estimation ($|E|$);
- Per-sample raw library size ($\sum_i K_{ij}$);
- Sample size factors ($s_j$);
- Per-sample normalized library size ($\sum_i \tilde{K}_{ij}$);
- Minimum and maximum size factors ($\min(s_j), \max(s_j)$);
- Max/min size-factor ratio ($\max(s_j) / \min(s_j)$);
- Relevant count and sparsity summaries (positive geometric mean count $|E^+|$, zero-fraction).

### 10.2 Configurable Threshold Policy

Threshold configurations in `NormalizationConfig` default to `None` (unset).
If an organization formally approves specific thresholds, breaches emit `configured_threshold_exceeded` (`Severity.REVIEW`), requiring a structured `NormalizationReviewDisposition`.

---

## 11. Data Structures and Schemas

```python
@dataclass(frozen=True, slots=True)
class SampleSizeFactor:
    sample_id: str
    size_factor: float
    raw_library_size: int
    normalized_library_size: float
    cohort_id: str
    response_eligibility: ResponseEligibility
    matched_control_sample_ids: tuple[str, ...]
    treatment_control_status: str

@dataclass(frozen=True, slots=True)
class NormalizedCohortResult:
    status: Status
    cohort: NormalizationCohort
    findings: tuple[Finding, ...]
    size_factors: tuple[SampleSizeFactor, ...]
    metrics: CohortNormalizationMetrics
    diagnostic_normalized_matrix: CountMatrix | None
    raw_count_matrix: CountMatrix
    method: str
    r_environment_info: dict[str, str]
    execution_provenance: dict[str, Any]
    review_disposition: NormalizationReviewDisposition | None = None
```

---

## 12. Verification Plan & Test Matrix

The implementation must include comprehensive tests covering:
1. **Normalization Cohort Construction:**
   - Treatment and matched controls preserved together despite differing `assay_id`;
   - Study and experiment boundary isolation (strictly no cross-study pooling);
   - Ambiguous technical context rejected with review.
2. **Raw-Count Immutability:**
   - Raw count values remain unchanged bit-for-bit.
3. **Official DESeq2 Execution & Numerical Validation:**
   - Execution against pinned R 4.4.3 / Bioconductor 3.20 / DESeq2 1.46.0;
   - Numerical equivalence against fixed reference fixtures;
   - Graceful detection and explicit error reporting if R is missing.
4. **Canonical Collisions:**
   - Unresolved collisions halt at `NEEDS_REVIEW` with `diagnostic_normalized_matrix = None` and raw matrix preserved.
   - Verification that no automatic collapsing (sum/mean/pick) takes place.
5. **Method Gating (Ratio vs. Poscounts):**
   - Standard ratio failure on all-zero geometric means emits `standard_size_factors_failed` and halts at `NEEDS_REVIEW`;
   - Poscounts executes only under explicit authorization.
6. **Missing Control Response Eligibility:**
   - Treatment samples without matched controls are tagged `ResponseEligibility.INELIGIBLE_MISSING_MATCHED_CONTROL`.
7. **Threshold Configuration:**
   - Unset (`None`) thresholds emit diagnostics without triggering review findings;
   - Configured threshold breaches emit structured review findings.
8. **Upstream QC/Harmonization Gating:**
   - Handled correctly for `PASS`, `PASS_WITH_WARNINGS`, authorized `NEEDS_REVIEW`, and `FAIL`.
9. **Full Regression Suite:**
   - Entire existing test suite (156 tests) must remain 100% green.
