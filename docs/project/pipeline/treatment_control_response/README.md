# Treatment-Control Response Construction

**Historical Identifier**: Step 7A
**Validation Status**: **VALIDATED** (Commit: `cd954f4`)

---

## 1. Scientific & Computational Purpose
Constructs scientifically valid, auditable treatment-versus-matched-control response records (`MolecularResponseContrast`). Evaluates biological replicate sufficiency, developmental stage compatibility, and organoid context matching without altering raw counts or leaking downstream DNT labels.

### Core Scientific Invariants:
1. **Explicit Matched Controls**: Treatment samples must be explicitly paired with matched vehicle controls from the same study and experiment.
2. **Developmental Stage Alignment**: Blocks comparisons across mismatched developmental ages (e.g. Day 6 vs Day 35).
3. **Biological Replicate Gating**: Requires $\ge 2$ true biological replicates per arm for inferential downstream differential expression. Technical replicates are audited but never inflate biological replicate counts.
4. **No Dose Pooling**: Distinct chemical concentrations or exposure durations produce separate contrast records.
5. **Zero Gene Filtering**: All genes remain intact for downstream testing.

---

## 2. Inputs & Outputs
* **Inputs**: `NormalizedDataset` (from Normalization) + canonical metadata mapping dictionary.
* **Outputs**: `ResponseContrastDataset` containing `MolecularResponseContrast` records, status (`ELIGIBLE`, `BLOCKED`, `NEEDS_REVIEW`), and replicate counts.

---

## 3. Code, Config & Test Map
* **Primary Code**:
  * [`src/response_builder/models.py`](../../../../src/response_builder/models.py): Contrast data models (`MolecularResponseContrast`, `ContrastExposure`).
  * [`src/response_builder/builder.py`](../../../../src/response_builder/builder.py): Contrast builder engine (`build_response_contrasts`).
* **Test Suite**:
  * [`tests/response_builder/test_builder.py`](../../../../tests/response_builder/test_builder.py)
  * [`tests/response_builder/test_response_models.py`](../../../../tests/response_builder/test_response_models.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: Contrast marked as `BLOCKED` with insufficient replicates**:
  * *First file to inspect*: [`src/response_builder/builder.py`](../../../../src/response_builder/builder.py) (inspect biological replicate counting).
  * *Second file to inspect*: Metadata sample definitions (verify `biological_replicate_id` distinguishes true biological replicates from technical replicates).
* **Symptom: Contrast blocked by developmental age mismatch**:
  * Verify `developmental_age` normalization in metadata (e.g. "Day 6" vs "d6").
