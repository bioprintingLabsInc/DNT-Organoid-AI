# Stage 7B — Validation Evidence Package & Benchmark Ledger

**Repository:** `bioprintingLabsInc/DNT-Organoid-AI`  
**Stage:** 7B (Bulk Treatment-versus-Matched-Control Differential Expression)  
**Validated Commit SHA:** `deb96dd27405c1ed85bbde78bda4704bf22e770c`  
**Validation Status:** **PASS** (100% agreement, zero numerical or design discrepancies)  
**Verification Date:** 2026-11-16  

---

## 1. Executive Summary

This validation evidence package provides permanent scientific and computational verification for **Stage 7B (Bulk Treatment-versus-Matched-Control Differential Expression)** of the DNT Organoid AI platform.

Stage 7B implements standard, metadata-governed differential gene expression analysis using Bioconductor DESeq2. To verify reproducibility and numerical accuracy, an independent direct R verification framework (`scripts/independent_deseq2_verification.py`) was developed and executed. This framework decouples from the application pipeline, constructs an isolated DESeq2 model in a standalone R process, and compares all statistical metrics gene-by-gene against the application pipeline output.

All acceptance criteria have passed without exception:
- **Zero Metadata Fabrication:** Complete elimination of synthetic biological replicates, default exposure parameters, and fallback design covariates (`source_1`, `batch_1`).
- **Mathematical Equivalence:** Numerical tolerance across all 6 statistical metrics strictly within $\le 5.0004 \times 10^{-9}$ (tolerance threshold: $1.0 \times 10^{-6}$).
- **Gene Universe Conservation:** 100% gene-for-gene concordance across all 60,207 canonical Ensembl GRCh38 release 112 genes.
- **Directional Concordance:** 100.0% sign agreement on treatment-vs-control fold changes (32,227 / 32,227 non-zero genes).
- **Statistical Design Verification:** Exact concordance in design formula (`~ condition`), design matrix rank (2), and residual degrees of freedom (6).

---

## 2. Software & Runtime Environment

The validation was executed and reproduced under the following pinned environment:

### Python Runtime
| Package / Runtime | Exact Version |
| :--- | :--- |
| **Python** | `3.14.6` |
| **pandas** | `3.0.3` |
| **numpy** | `2.5.1` |
| **scipy** | `1.18.0` |
| **pytest** | `9.1.1` |

### R & Bioconductor Runtime
| Software / Package | Exact Version | Binary Path |
| :--- | :--- | :--- |
| **R Engine** | `R version 4.4.3 (2025-02-28 ucrt)` | `C:\Users\chari\.r\R-4.4.3\bin\x64\Rscript.exe` |
| **Bioconductor** | `3.20` | Native |
| **DESeq2** | `1.46.0` | Bioconductor 3.20 release |

---

## 3. Test Suites & Execution Commands

### 3.1 Independent DESeq2 Verification CLI Suite (15 Tests)
Command:
```powershell
python -m pytest tests/differential_expression/test_independent_verification_cli.py
```
**Result: 15 passed in 438.72s (100% PASS)**

| Test # | Test Name | Target Invariant Verified | Result |
| :---: | :--- | :--- | :---: |
| 1 | `test_missing_biological_replicates_fails` | Rejects missing/fabricated biological replicate IDs | PASS |
| 2 | `test_missing_chemical_params_fails` | Rejects missing required exposure parameters | PASS |
| 3 | `test_simple_condition_design_runs` | Simple two-arm contrast (`~ condition`) matches direct R | PASS |
| 4 | `test_valid_paired_biological_source_design` | Paired biological source design (`~ biological_source + condition`) matches direct R | PASS |
| 5 | `test_valid_batch_adjusted_design` | Batch-adjusted design (`~ batch + condition`) matches direct R | PASS |
| 6 | `test_valid_combined_source_and_batch_design` | Combined source & batch design (`~ biological_source + batch + condition`) matches direct R | PASS |
| 7 | `test_confounded_statistical_design_fails` | Blocks confounded biological source/batch designs before execution | PASS |
| 8 | `test_numerical_vs_missing_value_mismatch_fails` | Detects missing-value mismatches between pipeline and direct R | PASS |
| 9 | `test_missing_gene_mismatch_fails` | Detects gene universe omission in direct R | PASS |
| 10 | `test_incorrect_statistical_design_fails` | Rejects design formula discrepancies between pipeline and R | PASS |
| 11 | `test_treatment_control_direction_mismatch_fails` | Catches inverted treatment-vs-control fold change directions | PASS |
| 12 | `test_custom_valid_dataset_passes` | Arbitrary valid RNA-seq count matrix & metadata runs cleanly | PASS |
| 13 | `test_missing_biological_source_identifier_rejected_when_required` | Rejects missing biological source ID when design includes `biological_source` | PASS |
| 14 | `test_missing_batch_identifier_rejected_when_required` | Rejects missing batch ID when design includes `batch` | PASS |
| 15 | `test_conflicting_metadata_identifiers_rejected` | Rejects conflicting identifiers across metadata alias columns | PASS |

### 3.2 Complete Differential Expression Suite (71 Tests)
Command:
```powershell
python -m pytest tests/differential_expression
```
**Result: 71 passed in 956.39s (100% PASS)**
- `test_design.py`: 18 passed
- `test_differential_expression.py`: 33 passed
- `test_independent_verification_cli.py`: 15 passed
- `test_models.py`: 5 passed

---

## 4. Benchmark Verification: GSE232218

### 4.1 Benchmark Metadata & Configuration
- **Accession:** NCBI GEO [GSE232218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232218)
- **Biological Context:** Day 6 Human Cerebral Cortex Organoids (SOSR-CO), cell line WA09
- **Perturbagen:** Valproic Acid (CID: 3121) at 200 µM nominal concentration in water vehicle
- **Exposure:** 48 hours continuous exposure
- **Sample Cohort:** 4 treated organoids (`Sample_3901-AT-12`, `Sample_3901-AT-17`, `Sample_3901-AT-3`, `Sample_3901-AT-7`) vs. 4 matched vehicle controls (`Sample_3901-AT-1`, `Sample_3901-AT-10`, `Sample_3901-AT-15`, `Sample_3901-AT-5`)
- **Benchmark Config File:** [`config/benchmarks/gse232218_vpa_200um.json`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/config/benchmarks/gse232218_vpa_200um.json)
- **Sample Metadata File:** [`config/benchmarks/gse232218_samples.csv`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/config/benchmarks/gse232218_samples.csv)

### 4.2 Benchmark Execution Command
```powershell
python scripts/independent_deseq2_verification.py --config config/benchmarks/gse232218_vpa_200um.json
```

### 4.3 Statistical Model & Design Concordance
| Property | Application Pipeline (Stage 7B) | Direct Independent R / DESeq2 | Concordance Status |
| :--- | :--- | :--- | :---: |
| **Model Formula** | `~ condition` | `~ condition` | **Exact Match** |
| **Model Matrix Rank** | 2 | 2 | **Exact Match** |
| **Residual Degrees of Freedom** | 6 | 6 | **Exact Match** |
| **Treatment Estimability** | Verified (`qr()$rank` increases) | Verified (`qr()$rank` increases) | **Exact Match** |

### 4.4 Gene Accounting & Universe Conservation
```
Total Canonical Ensembl GRCh38 Release 112 Genes:       60,207
├── Zero-Count Genes (baseMean == 0, p-values NA):        27,980 (Exact agreement)
└── Non-Zero Evaluated Genes (Wald Tested):              32,227 (Exact agreement)
    ├── Complete Numerical Match across 6 Metrics:        32,227 (100.0%)
    └── Direction Match (log2FC sign concordance):       32,227 / 32,227 (100.0%)
Discrepant Genes / Missing Value Mismatches:             0 (0.0%)
```

### 4.5 Numerical Verification Across 6 Statistical Metrics
Tolerance threshold set to $\varepsilon = 1.0 \times 10^{-6}$.

| Statistical Metric | Application Pipeline | Independent Direct R | Max Absolute Difference | Within Tolerance? |
| :--- | :--- | :--- | :---: | :---: |
| **Base Mean (`baseMean`)** | Direct computation | DESeq2 `results()` | $5.0004 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |
| **log2 Fold Change (`log2FoldChange`)** | Direct computation | DESeq2 `results()` | $4.9998 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |
| **LFC Standard Error (`lfcSE`)** | Direct computation | DESeq2 `results()` | $4.9998 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |
| **Wald Statistic (`stat`)** | Direct computation | DESeq2 `results()` | $5.0000 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |
| **p-value (`pvalue`)** | Direct computation | DESeq2 `results()` | $4.9995 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |
| **Adjusted p-value (`padj` / BH FDR)** | Direct computation | DESeq2 `results()` | $4.9948 \times 10^{-9}$ | **YES** ($< 10^{-6}$) |

### 4.6 Execution Timing
- Application Pipeline Execution: `451.17s`
- Standalone Independent R Execution: `26.34s`
- Total Benchmark Runtime: `477.51s`
- Overall Benchmark Status: **PASS**

---

## 5. Strict Scientific Invariants Enforced

1. **No Fabricated Biological Replicates:** Replicate IDs are never synthesized from sample identifiers. Every sample must declare genuine replicate metadata.
2. **No Assumed Experimental Metadata:** Default concentrations, exposure durations, developmental stages, and organoid contexts are prohibited. Missing parameters fail immediately.
3. **No Fabricated Design Covariates:** Covariates `biological_source` and `batch` are resolved strictly from user metadata. Silent fallback to `"source_1"` or `"batch_1"` has been completely removed.
4. **Estimability and Full Rank Enforcement:** Designs with complete confounding or linear dependence are rejected with `ContrastStatus.BLOCKED` before model fitting.
5. **Raw Integer Count Requirement:** DESeq2 requires raw, non-negative integer counts. Normalized floating-point values (TPM, FPKM, CPM) are rejected from DESeq2 and must be processed through non-count analysis routes.

---

## 6. Reproduction Instructions for External Researchers

To reproduce the independent numerical verification on any workstation:

1. **Clone the Repository and Checkout the Freeze Tag:**
   ```bash
   git clone https://github.com/bioprintingLabsInc/DNT-Organoid-AI.git
   cd DNT-Organoid-AI
   git checkout stage7b-de-v1.0-freeze
   ```

2. **Set Up Python Environment:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Verify R and DESeq2 Availability:**
   Ensure R $\ge 4.4$ and DESeq2 $\ge 1.46$ are installed:
   ```bash
   Rscript -e "library(DESeq2); cat('DESeq2 version:', as.character(packageVersion('DESeq2')), '\n')"
   ```

4. **Execute the Benchmark Verification:**
   ```bash
   python scripts/independent_deseq2_verification.py --config config/benchmarks/gse232218_vpa_200um.json
   ```

5. **Execute Verification on Custom User Data:**
   ```bash
   python scripts/independent_deseq2_verification.py \
     --counts path/to/raw_counts.csv \
     --samples path/to/sample_metadata.csv \
     --agent-name "CompoundX" \
     --concentration 50.0 \
     --developmental-age "Day 35" \
     --exposure-duration 24.0 \
     --output-dir verification_output/
   ```
