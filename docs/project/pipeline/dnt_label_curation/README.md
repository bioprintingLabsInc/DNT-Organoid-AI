# DNT Label Curation

**Historical Identifier**: Stage 9  
**Scientific Name**: Developmental Neurotoxicity Reference Label Curation  
**Status**: **STATUS: PLANNED — NO IMPLEMENTATION YET**

---

## 1. Scientific Objective
To assemble a multi-source, consensus developmental neurotoxicity ground-truth reference database (positive DNT toxicants, negative reference chemicals, and uncertain chemicals) based on regulatory guidelines (EPA, OECD, EFSA, NTP DNT-DIVER).

> [!CRITICAL]
> **Zero Leakage Rule**: Reference DNT labels are curated strictly independently of gene expression data. No model training or evaluation may infer labels from expression patterns.

---

## 2. Intended Inputs & Outputs
* **Intended Input**: Chemical regulatory dossiers, in vivo DNT test outcomes, and published literature evidence.
* **Intended Output**: Version-controlled chemical label catalog (`data/labels/`) supporting positive, negative, and unresolved/uncertain states.

---

## 3. Implementation Plan
* Implementation is planned under `src/labels/`.
* No code paths or algorithms are currently implemented.
