# Bulk RNA-seq Normalization

**Historical Identifier**: Step 6C  
**Validation Status**: **VALIDATED** (Commit: `5ce67cc`)  
**Umbrella Position**: Data Intake — Failure Domain 5

---

## 1. Scientific & Computational Purpose
Partitions harmonized datasets into metadata-defined normalization cohorts (by study, experiment, modality, assay) and calculates official Bioconductor **DESeq2 size factors** via the median-of-ratios method using an isolated headless R bridge. Computes diagnostic normalized counts for technical QC without altering raw integer counts intended for downstream differential expression.

---

## 2. Inputs & Outputs
* **Inputs**: `HarmonizedDataset` + validated metadata + configuration.
* **Outputs**: `NormalizedDataset` containing size factors (`SampleSizeFactor`), diagnostic normalized count matrix, and cohort metrics (`CohortNormalizationMetrics`).

---

## 3. Code, Config & Environment Map
* **Primary Code**:
  * [`src/normalization/cohort.py`](../../../../../src/normalization/cohort.py): Forms cohorts from metadata.
  * [`src/normalization/r_bridge.py`](../../../../../src/normalization/r_bridge.py): Subprocess bridge to Rscript.
  * [`src/normalization/normalizer.py`](../../../../../src/normalization/normalizer.py): Orchestration engine (`normalize_bulk_dataset`).
* **R Execution Engine**:
  * [`scripts/r/estimate_size_factors.R`](../../../../../scripts/r/estimate_size_factors.R): Headless R DESeq2 runner.
  * [`scripts/r/environment-r.yaml`](../../../../../scripts/r/environment-r.yaml), [`scripts/r/renv.lock`](../../../../../scripts/r/renv.lock): Environment lockfiles (R 4.4.3, DESeq2 1.46.0).
* **Configuration**:
  * [`config/bulk_normalization.yaml`](../../../../../config/bulk_normalization.yaml): Execution timeouts and threshold definitions.
* **Test Suite**:
  * [`tests/normalization/test_normalizer.py`](../../../../../tests/normalization/test_normalizer.py)
  * [`tests/normalization/test_r_bridge.py`](../../../../../tests/normalization/test_r_bridge.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: Size factor calculation fails or subprocess times out**:
  * *First file to inspect*: [`src/normalization/r_bridge.py`](../../../../../src/normalization/r_bridge.py).
  * *Second file to inspect*: [`scripts/r/estimate_size_factors.R`](../../../../../scripts/r/estimate_size_factors.R).
  * *Common cause*: Rscript executable is not on PATH, or all genes in a sample have zero counts.
