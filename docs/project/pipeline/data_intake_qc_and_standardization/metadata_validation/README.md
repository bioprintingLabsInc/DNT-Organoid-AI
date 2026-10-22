# Metadata Validation

**Historical Identifier**: Steps 3–4
**Validation Status**: **VALIDATED** (Commit: `455265c`)
**Umbrella Position**: Data Intake — Failure Domain 2

---

## 1. Scientific & Computational Purpose
Validates canonical experimental metadata records against [`config/metadata_schema.yaml`](../../../../../config/metadata_schema.yaml). Enforces structural validity, controlled vocabularies, foreign-key integrity between samples, assays, conditions, exposures, and experimental provenance without silently mutating or guessing missing values.

---

## 2. Inputs & Outputs
* **Inputs**: Metadata mapping dictionary conforming to canonical contract.
* **Outputs**: `ValidationResult` containing overall status (`VALID`, `VALID_WITH_WARNINGS`, `INVALID`) and deterministic `findings` list.

---

## 3. Code, Config & Test Map
* **Primary Code**:
  * [`src/metadata/schema.py`](../../../../../src/metadata/schema.py): YAML schema loader and indexer.
  * [`src/metadata/result.py`](../../../../../src/metadata/result.py): Finding, severity, and result dataclasses.
  * [`src/metadata/validator.py`](../../../../../src/metadata/validator.py): Validation engine.
* **Configuration**:
  * [`config/metadata_schema.yaml`](../../../../../config/metadata_schema.yaml): Canonical schema contract.
* **Documentation**:
  * [`docs/data_dictionary/METADATA_DATA_DICTIONARY.md`](../../../../data_dictionary/METADATA_DATA_DICTIONARY.md)
* **Test Suite**:
  * [`tests/metadata/test_validator.py`](../../../../../tests/metadata/test_validator.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: Metadata validation fails with `INVALID`**:
  * *First file to inspect*: [`src/metadata/validator.py`](../../../../../src/metadata/validator.py) (inspect `validate()`).
  * *Second file to inspect*: [`config/metadata_schema.yaml`](../../../../../config/metadata_schema.yaml) (inspect required fields).
  * *Common cause*: Missing `condition_id` on a sample record or invalid `biological_source_type` vocabulary.
