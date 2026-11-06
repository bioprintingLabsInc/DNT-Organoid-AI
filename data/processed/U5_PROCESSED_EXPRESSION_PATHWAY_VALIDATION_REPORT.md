# Phase U5 Targeted Scientific Validation Report
## Processed Continuous Bulk RNA-seq Pathway (U3.2, U3.3, U3.4) & Stage 8 Cross-Study Synthesis

**Date:** 2026-11-05  
**Repository:** `bioprintingLabsInc/DNT-Organoid-AI`  
**Branch:** `feature/stage8a2-cross-study-comparability`  
**Execution Environment:** Windows (x86_64), Python 3.14.6, R 4.4.3, Bioconductor limma 3.62.2  
**Validation Scope:** Targeted validation of the continuous processed expression pathway (Phases U3.2, U3.3, U3.4) on real GEO benchmark dataset **GSE187006**, and descriptive cross-study comparability (Stage 8A.2) against **GSE232218**.

> **Scope Limitation Notice:** This report validates solely the continuous processed expression ingestion and inference pathway. Raw FASTQ quantification, raw count intake (Step 6A), published differential expression table ingestion (U3.1), and single-cell/single-nucleus pseudobulk pipelines (Phase U4) remain separate milestones. **Phase U5 as an entirety is NOT declared complete.**

---

## 1. Executive Summary & Verification Matrix

| Requirement | Description | Status | Evidence / Artifact |
|---|---|:---:|---|
| **REQ-U5-1** | U3.4 format detection recognizes real TPM input and routes to U3.2 | **PASS** | `src/routing/detector.py`, `src/routing/dispatcher.py` |
| **REQ-U5-2** | U3.2 preserves floating-point continuous values, gene/sample identities, and blocks integer casting | **PASS** | `src/normalization/provenance_gate.py`, `src/ingestion/processed_matrix_reader.py` |
| **REQ-U5-3** | U3.3 executes official R/Bioconductor limma-trend for eligible biological replicates and blocks Python approximation | **PASS** | `src/differential_expression/limma_runner.py`, `scripts/r/run_limma_trend.R` |
| **REQ-U5-4** | Independent direct R reproduction with documented numerical discrepancies | **PASS** | `data/interim/direct_r_u1m_limma.csv`, max diff $\le 1.2 \times 10^{-8}$ on 1-to-1 loci |
| **REQ-U5-5** | Evidence files, verification suites, and SHA-256 output hashes generated | **PASS** | `data/processed/stage8a2_artifacts_sha256.json` |
| **REQ-U5-6** | Scientifically qualified Stage 8A.1 and Stage 8A.2 artifacts produced | **PASS** | Registries (U1M, U2F, Combined) & Cross-Study Comparison CSVs |
| **REQ-U5-7** | Preservation of GSE232218 results and explicit experimental regimen separation | **PASS** | Dose (200 µM vs 1 mM), stage (Day 6 vs Day 42), and duration (48h vs 72h) separated |
| **REQ-U5-8** | Strict scope bounding (processed pathway only; raw FASTQ, count, DEG, single-cell separate) | **PASS** | Explicitly stated and enforced |

---

## 2. Requirement Details & Empirical Evidence

### Requirement 1: U3.4 Format Detection & Route Dispatch
- **Objective:** Verify that U3.4 correctly recognizes real input files as processed continuous TPM matrices and routes to U3.2 without false classification as raw counts.
- **Observed Behavior:**
  - Asset: `data/raw/GSE187006/GSE187006_Processed_data.csv.gz` (17,882 rows × 12 columns, $\log_2(\text{TPM})$ values).
  - `detector.detect_format()` established:
    - `detected_format`: `processed_expression_matrix`
    - `confidence`: `Confidence.HIGH`
    - `evidence`: `('gzip-compressed tabular structure', '50 inspected data rows', 'non-integer numeric values')`
    - `findings`: `()` (zero warnings, zero errors)
  - `router._select_route()` returned `Route.PROCESSED_EXPRESSION_PIPELINE`.
  - `dispatcher.plan_route_execution()` generated `RouteExecutionPlan(route=Route.PROCESSED_EXPRESSION_PIPELINE, status="READY_FOR_EXECUTION", executor_id="processed_expression_pipeline_runner", target_format="differential_expression_dataset")`.
- **Verdict:** **PASS**

---

### Requirement 2: U3.2 Intake & Provenance Gating
- **Objective:** Verify that U3.2 preserves continuous TPM values, gene identities, and sample identifiers without treating continuous values as raw integer counts.
- **Observed Behavior:**
  - Matrix parsed: 17,882 rows (gene symbols) × 12 samples (`U1M.NC1..3`, `U1M.V1..3`, `U2F.NC2..4`, `U2F.V2..4`).
  - Gene Harmonization:
    - 15,572 genes uniquely mapped to canonical Ensembl Release 112 IDs.
    - 2,310 genes excluded (unmapped or colliding aliases).
    - **Law of Conservation of the Gene Universe:** $15,572 + 2,310 = 17,882$ exactly.
  - Scale Validation: Provenance established as `ExpressionUnit.TPM`, `TransformationScale.LOG2`, `is_log_transformed=True`. Negative values (present for $0 < \text{TPM} < 1.0$) accepted as valid fractional TPM logs.
  - Scientific Isolation: DESeq2 intake was completely blocked (`deseq2_eligible=False`, `blocked_analysis_engines=('DESeq2',)`). Zero rounding or integer conversion was applied.
  - Intake Quality Control: `validate_processed_matrix` returned status `PASS` with 0 QC findings.
- **Verdict:** **PASS**

---

### Requirement 3: U3.3 limma-trend Inference via Official R/Bioconductor
- **Objective:** Verify that U3.3 invokes the official Bioconductor `limma-trend` empirical Bayes pipeline via Rscript, verifies biological replication ($n \ge 2$), and blocks Python approximation from entering accepted results.
- **Observed Behavior:**
  - R Environment: Official R 4.4.3 with Bioconductor `limma` v3.62.2.
  - Replicate Gating:
    - Contrast 1 (`U1M_VPA_1mM_72h_vs_NC`): $n=3$ treated vs $n=3$ vehicle control ($n \ge 2$, verified independent organoid replicates) $\to$ Status `ELIGIBLE`.
    - Contrast 2 (`U2F_VPA_1mM_72h_vs_NC`): $n=3$ treated vs $n=3$ vehicle control ($n \ge 2$, verified independent organoid replicates) $\to$ Status `ELIGIBLE`.
  - Python Approximation: When R is missing, `_execute_python_fallback` flags `result_status="UNSUPPORTED_APPROXIMATION"` and contrast status `BLOCKED`. It is strictly prohibited from entering Stage 8A registries.
  - Real Execution: Both contrasts were executed using `run_limma_trend.R`:
    - Model formula: `~ 0 + condition` (rank = 2, residual degrees of freedom = 4).
    - Empirical Bayes moderation: `eBayes(..., trend=TRUE)`.
    - P-value adjustment: Benjamini–Hochberg (BH).
    - Moderated t-statistic preserved, Wald standard error uninvented (`None`).
    - Results:
      - U1M: 15,561 evaluated canonical genes; 3,540 with FDR < 0.05. Top hit: `HLA-B` ($\log_2\text{FC}=2.858$, $t=27.44$, $p=8.10 \times 10^{-9}$, $p_{\text{adj}}=1.26 \times 10^{-4}$).
      - U2F: 15,561 evaluated canonical genes; 3,118 with FDR < 0.05. Top hits: `PRKCB` ($\log_2\text{FC}=1.891$, $t=18.67$, $p=2.00 \times 10^{-8}$), `HLA-B` ($\log_2\text{FC}=3.069$, $t=16.52$, $p=5.74 \times 10^{-8}$).
- **Verdict:** **PASS**

---

### Requirement 4: Independent Numerical Reproduction & Discrepancy Audit
- **Objective:** Independently reproduce differential expression using direct R execution on the same inputs, compare against pipeline outputs, and document numerical discrepancies.
- **Comparison Method:**
  - Direct independent R script executed `limma::lmFit`, `limma::makeContrasts`, `limma::contrasts.fit`, `limma::eBayes(trend=TRUE)`, `limma::topTable` on the raw expression CSV.
  - Compared 15,004 shared gene symbols between direct R and pipeline results.
- **Discrepancy Findings:**
  - Mean absolute difference in $\log_2\text{FC}$: $1.88 \times 10^{-5}$
  - Mean absolute difference in moderated $t$: $2.23 \times 10^{-2}$
  - For all 1-to-1 non-colliding loci: max absolute $\log_2\text{FC}$ difference is $\le 1.20 \times 10^{-8}$ (exact to the 8-decimal JSON serialization limit).
  - **Identified Discrepancy Root Cause:** A single symbol (`MIR222HG`) showed a $\log_2\text{FC}$ discrepancy of 0.282. Investigation revealed that in the source data, both `LINC02595` (row 9251) and `MIR222HG` (row 9978) were present and both map to Ensembl ID `ENSG00000270069`. The pipeline strictly enforced canonical collision deduplication (retaining the first row and excluding duplicate entries to preserve the canonical gene universe), whereas direct unharmonized R evaluated both symbol rows independently.
- **Verdict:** **PASS**

---

### Requirement 5: Evidence Artifacts & Cryptographic Checksums
- **Objective:** Persist all generated artifacts and record full-file SHA-256 checksums.
- **Artifact Checksums (`stage8a2_artifacts_sha256.json`):**
  - `GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv`: `8016eb82af04e0821a34f56a634151dfd5b5228813c24659450e4365c72109ab`
  - `GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv`: `0128eecff99a1715d6ae2bfdd851ec8fb9460cd2d6ca77ae588190a72558a3a1`
  - `GSE187006_VPA_1mM_condition_gene_response_registry.csv`: `37ad514ada52832adf240509345b14b094821a05064599403812f0f40b4cca7c`
  - `GSE232218_vs_GSE187006_U1M_cross_study_comparison.csv`: `e51677d234da0062e934bd3ba57e7a995221e3c2f2cb319522d1487f7d87d988`
  - `GSE232218_vs_GSE187006_U2F_cross_study_comparison.csv`: `4cc2b8cb304eaee6c40d29b4adcd99e395de5b65872ae9d257dd307821fc9c41`
  - Source data `GSE187006_Processed_data.csv.gz`: `948795464348bc692ddfe97d617e4f2261616e13fb1d19abbd3d733b3faf6e00`
  - Baseline `GSE232218_VPA_200uM_condition_gene_response_registry.csv`: `67efb598d1a1dd15c2d12f6e9bc3537fb117b9b66a5e12ec6850c978b77a7f6f`
- **Verdict:** **PASS**

---

### Requirement 6 & 7: Stage 8A.1 & 8A.2 Scientific Qualification & Context Preservation
- **Objective:** Generate governed condition-level registries and perform descriptive cross-study comparison against GSE232218 while strictly preserving experimental differences and prohibiting invalid meta-analysis pooling.
- **Stage 8A.1 Registries:**
  - `GSE187006_U1M_VPA_1mM`: 15,561 records, 0 numerical mismatches, 0 leakage violations.
  - `GSE187006_U2F_VPA_1mM`: 15,561 records, 0 numerical mismatches, 0 leakage violations.
  - Combined `GSE187006_VPA_1mM`: 31,122 records, 0 numerical mismatches, 0 leakage violations.
- **Stage 8A.2 Comparability Assessment:**
  - Evaluated pairs:
    1. GSE232218 (200 µM, 48h, Day 6 SOSR-CO WA09) vs GSE187006 U1M (1 mM, 72h, Day 42 forebrain organoids)
    2. GSE232218 (200 µM, 48h, Day 6 SOSR-CO WA09) vs GSE187006 U2F (1 mM, 72h, Day 42 forebrain organoids)
  - Comparability Gate Results:
    - Compound Match: Valproic Acid / CID:3121 (**MATCH**)
    - Dose: 200 µM vs 1.0 mM (**MISMATCH**, 5-fold concentration disparity)
    - Developmental Stage: Day 6 neural induction vs Day 42 mature organoids (**MISMATCH**)
    - Exposure Duration: 48h vs 72h (**MISMATCH**)
    - Organoid Platform: Cerebral Cortex Organoid (SOSR-CO) vs Forebrain Organoid (**MISMATCH**)
    - Biological Replication: Both studies verified ($n=4$ vs $n=3$, fully replicated)
    - **Synthesis Decision:** `is_synthesis_justified = False`, `cross_study_status = CrossStudyStatus.NONCOMPARABLE_CONTEXT`.
    - **Synthesis Boundary Enforcement:** Across all 60,221 evaluated genes, `pooled_log2_fold_change` is strictly `None` and `pooled_standard_error` is strictly `None`. Unweighted averaging was completely blocked.
- **Descriptive Concordance Metrics:**
  - Total union gene universe: 60,221 loci
  - Shared canonical genes with valid estimates: 15,547 loci
  - **GSE232218 vs GSE187006 U1M:**
    - Concordant UP: 4,877 genes
    - Concordant DOWN: 4,739 genes
    - Discordant: 5,736 genes
    - Undetermined: 44,869 genes (unexpressed or absent in one study)
    - Concordance Rate (directional genes): **62.64%** (9,616 / 15,352)
  - **GSE232218 vs GSE187006 U2F:**
    - Concordant UP: 4,632 genes
    - Concordant DOWN: 4,862 genes
    - Discordant: 5,858 genes
    - Undetermined: 44,869 genes
    - Concordance Rate (directional genes): **61.84%** (9,494 / 15,352)
- **Verdict:** **PASS**

---

### Requirement 8: Explicit Scope Bounding
- **Scope Boundary Statement:**
  - This validation strictly pertains to the continuous processed expression matrix intake, provenance gating, and limma-trend differential inference pathway.
  - Raw FASTQ alignment (Phase U1/F1.1), raw count matrix ingestion (Step 6A), published DEG summary tables (U3.1), and single-cell/single-nucleus pseudobulk inference (Phase U4) are maintained as independent pathways.
  - Universal Ingestion Phase U5 is NOT declared complete; only the continuous processed expression subsystem is validated.
- **Verdict:** **PASS**
