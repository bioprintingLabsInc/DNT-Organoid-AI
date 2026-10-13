# Compound-Gene & Pathway Relationships

**Historical Identifier**: Stage 8D  
**Scientific Name**: Compound-Gene and Compound-Pathway Relationships  
**Status**: **STATUS: PLANNED — NO IMPLEMENTATION YET**

---

## 1. Scientific Objective & Primary Purpose

The primary purpose of Compound-Gene & Pathway Relationships is to synthesize experimental compound-perturbation associations directly observed in human brain organoids:
* **Differential Molecular Response** $\to$ Experimental **compound–gene relationships**
* **Pathway & Signaling Features** $\to$ Experimental **compound–pathway relationships**

### Required Context Attributes to Preserve:
For every observed compound–gene and compound–pathway relationship, the representation must preserve:
* **Compound Identity**: Chemical name, CAS number, PubChem CID;
* **Target Entity**: Canonical gene ID / symbol or pathway identifier / name;
* **Perturbation Metrics**: Effect magnitude ($\log_2\text{FC}$ or GSEA NES), response direction (up/down/positive/negative), statistical confidence / uncertainty ($p$-value, FDR, SE);
* **Experimental Regimen**: Concentration / dose, exposure duration, organoid developmental stage;
* **Study & Context**: Study ID, experiment ID, biological source / organoid line, reproducibility metrics across replicates;
* **Provenance**: Upstream execution hashes and contrast IDs.

---

## 2. External Toxicological Knowledge & Circularity Safeguards

External chemical databases (e.g., ChEMBL, PubChem, ToxCast, AOP-Wiki, or known literature mechanisms-of-action):
* **Must be treated strictly as a separate annotation and post-hoc validation layer**;
* **Must NOT be automatically incorporated as v1 predictive input features**.

> [!WARNING]
> **Circularity & Leakage Risk**:
> Many external chemical-pathway databases incorporate historical neurotoxicity annotations, assay readouts, or expert consensus that overlap directly with the regulatory DNT endpoints our models aim to predict. Blending external toxicology labels into feature vectors creates severe **circularity**, artificially inflating model performance and masking whether human brain organoids actually generate independent predictive transcriptomic signal.

---

## 3. Implementation Plan

* Implementation will reside in `src/features/compound_pathway/` and `tests/features/compound_pathway/`.
* No code paths or algorithms are currently implemented.
