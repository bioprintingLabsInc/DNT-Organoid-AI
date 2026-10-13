# Git Branch, Tag, and Release Map

This document explains the repository branching strategy, validated release milestones, frozen tags, and future branch naming conventions.

---

## 1. Active Branches and Linage

| Branch Name | Role & Scope | HEAD Commit SHA | Validation Status |
| :--- | :--- | :--- | :--- |
| **`main`** | Production integration trunk | `f1a7ad8` | Contains validated Stages 1–7B (Differential Expression v1) |
| **`feature/streamlit-demo-hardening`** | Web interface hardening & GEO resolution | `ca26d34` | Validated Streamlit web interface on GSE232218 |
| **`feature/stage8c-pathway-analysis-v1`** | Pathway and signaling features v1 | `7a69c25` | **FROZEN RELEASE BASELINE** (`stage8c-v1-frozen`) |
| **`chore/repository-organization-handoff-v1`** | Documentation & handoff reorganization | Current | Documentation-only branch based on `stage8c-v1-frozen` |

### Historical Feature Branches (Archived for Traceability):
* `feature/gene-harmonization-v1` (commit `b904a41`): Step 6B Gene Harmonization locked.
* `feature/bulk-normalization-v1` (commit `5ce67cc`): Step 6C Bulk Normalization locked.
* `feature/contrast-builder-v1` (commit `cd954f4`): Step 7A Treatment-Control Contrast Builder locked.
* `feature/bulk-differential-expression-v1` (commit `f1a7ad8`): Step 7B Bulk Differential Expression locked.
* `feature/streamlit-analysis-app` (commit `8ab1503`): Initial web interface implementation.

---

## 2. Immutable Release Tags

### Release Tag: `stage8c-v1-frozen`
* **Commit SHA**: `7a69c259b9be1ea9cff55dd76b8db365df105b9f`
* **Tag Type**: Annotated Git Tag (`refs/tags/stage8c-v1-frozen`)
* **Scientific Milestone**: Pathway and Signaling Feature Foundation (GSEA GO BP and KEGG) verified on GSE232218.
* **Freeze Invariant**: This commit is immutable. No downstream feature work or refactoring may alter this tag or its referenced tree.

---

## 3. Future Branch Naming Standards

To maintain scientific clarity, future branches must use **scientific functional terminology** rather than historical stage numbers:

| Planned Scientific Milestone | Preferred Future Branch Name | Deprecated Historical Branch Name |
| :--- | :--- | :--- |
| **Cross-Study Gene Response Features** | `feature/cross-study-gene-response` | `feature/stage8a` |
| **Gene Coresponse & Network Features** | `feature/gene-coresponse-networks` | `feature/stage8b` |
| **Compound-Gene-Pathway Relationships** | `feature/compound-pathway-relationships` | `feature/stage8d` |
| **DNT Label Curation** | `feature/dnt-label-curation` | `feature/stage9` |
| **Compound-Safe Dataset Splitting** | `feature/compound-safe-splitting` | `feature/stage10` |
| **DNT Signature Discovery** | `feature/dnt-signature-discovery` | `feature/stage11` |
| **Model Training & Calibration** | `feature/model-training-calibration` | `feature/stage12` |

---

## 4. Release Governance Rules

1. **Tag Immobility**: Once an annotated release tag (e.g. `stage8c-v1-frozen`) is published to `origin`, it must never be deleted, moved, or recreated.
2. **Acceptance Before Merging**: A branch must pass 100% of repository tests (`pytest -q`), pass `git diff --check`, and reproduce benchmark outputs before merging into `main`.
3. **No Unrecorded Code Modifications**: Any refactoring must be accompanied by explicit test updates and documented in `docs/project/RELEASE_HISTORY.md`.
