# Pathway & Signaling Features

**Historical Identifier**: Stage 8C
**Scientific Name**: Pathway and Signaling Feature Foundation
**Release Status**: **FROZEN v1**
**Immutable Release Tag**: [`stage8c-v1-frozen`](https://github.com/bioprintingLabsInc/DNT-Organoid-AI/tree/stage8c-v1-frozen)
**Frozen Baseline Commit**: `7a69c259b9be1ea9cff55dd76b8db365df105b9f`

---

## 1. Scientific & Computational Purpose
Converts the unbiased genome-wide treatment-versus-matched-control molecular response from bulk RNA-seq (produced in Stage 7B) into standardized, interpretable pathway- and signaling-level feature representations using preranked Gene Set Enrichment Analysis (GSEA) across **Gene Ontology Biological Process (GO BP)** and **KEGG Pathways**.

> [!IMPORTANT]
> **Foundational Governance Principle**: This module operates **without using any developmental neurotoxicity (DNT) outcome labels or expected biological answers**. It strictly characterizes the cellular perturbation phenotype of the chemical exposure in human brain organoids.

---

## 2. Frozen Gene Universe & Reconciliation Accounting

Audited directly on the public **GSE232218** Valproic Acid 200 µM pilot dataset:

```
Stage 7 Canonical Evaluated Genes                    60,207
├── All-Zero Count Genes (baseMean == 0)            -27,980
└── Genes with Statistically Tested Wald             32,227
    ├── Cook's Distance Outliers                        -87
    └── STATISTICALLY ELIGIBLE UNIVERSE              32,140 (100% Preserved)
        ├── Entrez-Unmapped Valid Genes              -7,122 (Preserved in non-Entrez universe)
        └── Entrez-Mapped Valid Genes                25,018
            └── Duplicate Entrez Resolved                -2 (MAX_ABS_WALD resolution)
                └── FINAL GSEA RANKED GENES          25,016
```

### Architectural Preservation of the 7,122 Non-Entrez Genes:
\[
\boxed{\text{No statistically valid gene is discarded just because it lacks an Entrez ID}}
\]
The 7,122 statistically eligible genes lacking NCBI Entrez IDs in Ensembl 112 are **not QC failures**. They are fully preserved in [`data/processed/GSE232218_VPA_200uM_retained_non_entrez_genes.csv`](../../../../data/processed/GSE232218_VPA_200uM_retained_non_entrez_genes.csv) with `future_feature_eligibility=True` and remain available for downstream cross-study gene-response modeling.

### Exact Biological Duplicate Entrez Resolution (`MAX_ABS_WALD`):
* **Target 1: Entrez ID `125775236`**:
  * Candidate 1: `ENSG00000240758` (`HILPDA-AS1`, lncRNA), Wald = `+1.258533`
  * Candidate 2: `ENSG00000272601` (`EFCAB3P1`, pseudogene), Wald = `-1.133894`
  * **Selected Representative**: `ENSG00000240758` (`HILPDA-AS1`) because $|1.258533| > |-1.133894|$.
* **Target 2: Entrez ID `100874204`**:
  * Candidate 1: `ENSG00000250012` (`SLC41A3-AS1`, lncRNA), Wald = `-0.387471`
  * Candidate 2: `ENSG00000250218` (`ALDH1L1-AS1`, lncRNA), Wald = `+0.240024`
  * **Selected Representative**: `ENSG00000250012` (`SLC41A3-AS1`) because $|-0.387471| > |0.240024|$.
* *Note on Test Fixtures*: Mentions of `7157` (`TP53`) and `101` (`GENE1`) in earlier notes were synthetic mock test fixtures only, not real pilot duplicates.

---

## 3. Pathway Database Results

* **GO Biological Process**: 4,893 gene sets evaluated ($15 \le \text{size} \le 500$); 1,151 significant at BH-adjusted $p < 0.05$ (23.52%).
* **KEGG Pathways**: 348 pathways evaluated ($15 \le \text{size} \le 500$); 129 significant at BH-adjusted $p < 0.05$ (37.07%).

---

## 4. Code, Reference & Test Map

* **Primary Code**:
  * [`src/features/pathway/mapping.py`](../../../../src/features/pathway/mapping.py): Maps Ensembl canonical IDs to Entrez IDs.
  * [`src/features/pathway/ranking.py`](../../../../src/features/pathway/ranking.py): Builds ranked gene set and resolves duplicate Entrez targets.
  * [`src/features/pathway/models.py`](../../../../src/features/pathway/models.py): Feature, ranked set, and enrichment models.
  * [`src/features/pathway/qc.py`](../../../../src/features/pathway/qc.py): Mathematical QC reconciliation and metrics.
  * [`src/features/pathway/provenance.py`](../../../../src/features/pathway/provenance.py): Assembles full execution provenance audit.
  * [`src/features/pathway/r_bridge.py`](../../../../src/features/pathway/r_bridge.py): Headless R bridge executing clusterProfiler 4.14.6.
  * [`src/features/pathway/pipeline.py`](../../../../src/features/pathway/pipeline.py): End-to-end orchestration pipeline.
* **References & Scripts**:
  * [`data/references/kegg_hsa_pathways.rds`](../../../../data/references/kegg_hsa_pathways.rds): Offline KEGG snapshot.
  * [`scripts/r/run_preranked_gsea.R`](../../../../scripts/r/run_preranked_gsea.R): Headless R GSEA script.
  * [`scripts/run_vpa_pilot_stage8c.py`](../../../../scripts/run_vpa_pilot_stage8c.py): Verification execution script.
* **Test Suite**:
  * [`tests/features/pathway/test_statistically_eligible_universe.py`](../../../../tests/features/pathway/test_statistically_eligible_universe.py)
  * [`tests/features/pathway/test_ranking.py`](../../../../tests/features/pathway/test_ranking.py)
  * [`tests/features/pathway/test_mapping.py`](../../../../tests/features/pathway/test_mapping.py)
  * [`tests/features/pathway/test_gsea_integration.py`](../../../../tests/features/pathway/test_gsea_integration.py)
  * [`tests/features/pathway/test_qc_provenance.py`](../../../../tests/features/pathway/test_qc_provenance.py)
  * [`tests/features/pathway/test_leakage.py`](../../../../tests/features/pathway/test_leakage.py)
