# Practical Platform Troubleshooting Guide

This guide enables developers, scientists, and AI agents to quickly diagnose and resolve pipeline failures without prior project context.

---

## 1. Quick Diagnostic Matrix

| Failure Symptom | Likely Failure Domain | First File to Inspect | Second File to Inspect | Relevant Config / References | Relevant Test Suite | Escalation Path |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **"Unrecognized input data format" or Empty Matrix Error** | Intake Domain 1: Ingestion & Format Detection | `src/routing/detector.py` | `src/ingestion/inspection.py` | `src/ingestion/models.py` | `tests/routing/test_routing.py` | Check file magic bytes; ensure UTF-8 encoding and valid delimiters. |
| **"Canonical metadata is invalid; routing cannot be approved"** | Intake Domain 2: Metadata Validation | `src/metadata/validator.py` | `src/metadata/schema.py` | `config/metadata_schema.yaml` | `tests/metadata/test_validator.py` | Inspect `_findings` array in `ValidationResult`; check missing foreign keys. |
| **"Matrix sample IDs do not match metadata" or Low Spearman Correlation** | Intake Domain 3: Raw-Count QC | `src/qc/matrix.py` | `src/qc/replicates.py` | `config/qc_rules.yaml` | `tests/qc/test_validator.py` | Check sample header alignment; verify integer counts; inspect replicate IDs. |
| **High Unmapped Gene Fraction (>20%) or Collision Error** | Intake Domain 4: Gene Harmonization | `src/gene_harmonization/type_detector.py` | `src/gene_harmonization/harmonizer.py` | `config/gene_harmonization.yaml`<br>`data/references/ensembl_human_genes_v112.tsv` | `tests/gene_harmonization/test_harmonizer.py` | Verify whether input IDs are Ensembl, HGNC, or Entrez; check version suffix stripping (`.13`). |
| **Size Factors All Zero, Negative, or Subprocess Timeout** | Intake Domain 5: Normalization | `src/normalization/cohort.py` | `src/normalization/r_bridge.py` | `config/bulk_normalization.yaml`<br>`scripts/r/estimate_size_factors.R` | `tests/normalization/test_normalizer.py` | Ensure samples have non-zero shared counts; verify `Rscript` binary is on PATH; check R package versions (`DESeq2 1.46.0`). |
| **"insufficient_biological_replicates" or Age Mismatch Error** | Molecular Response: Contrast Construction | `src/response_builder/builder.py` | `src/response_builder/models.py` | `config/metadata_schema.yaml` | `tests/response_builder/test_builder.py` | Ensure at least 2 true biological replicates per arm; check developmental age normalization (e.g. "Day 6"). |
| **"Gene universe not conserved" or Model Fit Error** | Molecular Response: Differential Expression | `src/differential_expression/pipeline.py` | `src/differential_expression/r_bridge.py` | `scripts/r/run_deseq2_contrast.R` | `tests/differential_expression/test_differential_expression.py` | Check that evaluated + excluded genes equals total input rows; inspect R glm convergence warnings. |
| **"Offline KEGG pathways RDS not found" or clusterProfiler Timeout** | Features: Pathway & Signaling Analysis | `src/features/pathway/r_bridge.py` | `src/features/pathway/ranking.py` | `data/references/kegg_hsa_pathways.rds`<br>`scripts/r/run_preranked_gsea.R` | `tests/features/pathway/test_gsea_integration.py` | Confirm `kegg_hsa_pathways.rds` exists; check R memory limit; inspect tie-breaking logic in `GeneRankingBuilder`. |

---

## 2. Deep-Dive Runbooks by Failure Domain

### Domain 1: File Ingestion & Format Detection
* **Trigger**: A raw count table is uploaded, but status is `INVALID` or `MANUAL_REVIEW`.
* **Inspection Protocol**:
  1. Open `src/routing/detector.py`: check `detect_format()`.
  2. Inspect the first 5 lines of the data file. Does it have a header? Is the first column gene IDs?
  3. Ensure count values are positive integers without decimal points (e.g. `12` not `12.0`).

### Domain 2: Metadata Schema Validation
* **Trigger**: `MetadataValidator().validate(metadata)` produces errors.
* **Inspection Protocol**:
  1. Print findings: `[print(f.path, f.message) for f in result.findings if f.severity == 'ERROR']`.
  2. Common causes: missing `condition_id`, duplicate `sample_id`, invalid vocabulary term in `biological_source_type`.
  3. Compare with `config/metadata_schema.yaml`.

### Domain 3: Quality Control & Sample Reconcile
* **Trigger**: `validate_bulk_counts()` fails or issues warnings.
* **Inspection Protocol**:
  1. Check `DatasetResult.findings`.
  2. Are sample names in the matrix identical to `sample_id` in metadata?
  3. Are there samples with zero total reads?

### Domain 4: Gene Harmonization
* **Trigger**: `HarmonizedDataset.status` is `NEEDS_REVIEW` or `INVALID`.
* **Inspection Protocol**:
  1. Open `src/gene_harmonization/type_detector.py`. Run a small test lookup on 5 problematic IDs.
  2. Check if the gene is novel or non-canonical (unsupported in Ensembl 112).
  3. Ensure unmapped genes are preserved with status `UNMAPPED` rather than crashed.

### Domain 5: Headless R Execution Bridges
* **Trigger**: `REnvironmentInfo.is_available == False` or subprocess timeout.
* **Inspection Protocol**:
  1. Check `which Rscript` (or `where Rscript` on Windows).
  2. Run `Rscript -e "packageVersion('DESeq2'); packageVersion('clusterProfiler')"` to verify versions.
  3. If packages are missing, restore the R environment: `Rscript -e "renv::restore()"` using `scripts/r/renv.lock`.
