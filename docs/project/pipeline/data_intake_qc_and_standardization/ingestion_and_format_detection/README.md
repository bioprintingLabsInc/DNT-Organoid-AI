# Ingestion & Format Detection

**Historical Identifier**: Step 5
**Validation Status**: **VALIDATED** (Commit: `1bf8ddb`)
**Umbrella Position**: Data Intake — Failure Domain 1

---

## 1. Scientific & Computational Purpose
Determines the physical format and RNA-seq modality of incoming data assets (FASTQ, raw counts, UMI counts, H5AD, RDS) without altering raw content. Routes assets strictly to modality-appropriate pipelines (`BULK_RAW_COUNT_PIPELINE`, etc.) and flags declared/detected mismatches for review.

---

## 2. Inputs & Outputs
* **Inputs**: File path on filesystem (`.csv`, `.tsv`, `.txt`, `.fastq`, `.fastq.gz`, `.h5ad`, `.rds`) + declared metadata.
* **Outputs**: `AssetResult` tuple containing `detected_format`, `confidence`, `route`, `status` (`READY`, `NEEDS_REVIEW`, `INVALID`), and `findings`.

---

## 3. Code, Config & Test Map
* **Primary Code**:
  * [`src/ingestion/models.py`](../../../../../src/ingestion/models.py): Asset status enums, detection containers.
  * [`src/ingestion/inspection.py`](../../../../../src/ingestion/inspection.py): Header and magic byte inspection.
  * [`src/ingestion/registry.py`](../../../../../src/ingestion/registry.py): Registration against metadata.
  * [`src/routing/detector.py`](../../../../../src/routing/detector.py): Format detector rules.
  * [`src/routing/router.py`](../../../../../src/routing/router.py): Modality routing engine.
* **Test Suite**:
  * [`tests/ingestion/test_registration.py`](../../../../../tests/ingestion/test_registration.py)
  * [`tests/routing/test_routing.py`](../../../../../tests/routing/test_routing.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: File cannot be read or detected as `INVALID`**:
  * *First file to inspect*: [`src/routing/detector.py`](../../../../../src/routing/detector.py) (inspect `detect_format()`).
  * *Second file to inspect*: [`src/ingestion/inspection.py`](../../../../../src/ingestion/inspection.py) (inspect byte reader).
* **Symptom: Route assigned to `MANUAL_REVIEW`**:
  * Check for declared vs detected mismatch (e.g. declared `raw_counts` but detected floating-point expression).
