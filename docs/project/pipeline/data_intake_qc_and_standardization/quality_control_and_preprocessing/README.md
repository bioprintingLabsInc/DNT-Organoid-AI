# Quality Control & Preprocessing

**Historical Identifier**: Step 6A
**Validation Status**: **VALIDATED** (Commit: `3a8a7ef`)
**Umbrella Position**: Data Intake — Failure Domain 3

---

## 1. Scientific & Computational Purpose
Evaluates bulk RNA-seq count matrix structural validity and technical characteristics prior to gene harmonization and normalization. Calculates sample library sizes, zero-count fractions, and biological replicate Spearman correlations. Flags anomalous samples without silently discarding them.

---

## 2. Inputs & Outputs
* **Inputs**: Path to raw-count matrix (`.csv`, `.tsv`, `.txt`) + validated canonical metadata dictionary.
* **Outputs**: `DatasetResult` containing sample metrics (`SampleResult`), condition replicate correlations (`ConditionResult`), and rule findings.

---

## 3. Code, Config & Test Map
* **Primary Code**:
  * [`src/qc/matrix.py`](../../../../../src/qc/matrix.py): Count matrix reader and structural validation.
  * [`src/qc/metadata.py`](../../../../../src/qc/metadata.py): Reconciles matrix columns with metadata samples.
  * [`src/qc/metrics.py`](../../../../../src/qc/metrics.py): Library size, detected genes, zero fractions.
  * [`src/qc/replicates.py`](../../../../../src/qc/replicates.py): Average-rank Spearman replicate correlation.
  * [`src/qc/rules.py`](../../../../../src/qc/rules.py): Structural rules and warning threshold evaluations.
  * [`src/qc/validator.py`](../../../../../src/qc/validator.py): Orchestration engine (`validate_bulk_counts`).
* **Configuration**:
  * [`config/qc_rules.yaml`](../../../../../config/qc_rules.yaml): QC structural rules and thresholds.
* **Test Suite**:
  * [`tests/qc/test_matrix.py`](../../../../../tests/qc/test_matrix.py), [`tests/qc/test_validator.py`](../../../../../tests/qc/test_validator.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: QC status `FAIL` or `REVIEW`**:
  * *First file to inspect*: [`src/qc/validator.py`](../../../../../src/qc/validator.py).
  * *Second file to inspect*: [`src/qc/matrix.py`](../../../../../src/qc/matrix.py) (check non-integer or negative count errors).
  * *Common cause*: Matrix sample names do not match metadata IDs, or counts contain non-numeric tokens.
