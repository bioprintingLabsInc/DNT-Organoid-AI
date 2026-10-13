# Platform Release History & Historical Stage Crosswalk

This document records the chronological development milestones of the DNT Organoid AI platform and maps legacy engineering stage numbers to canonical scientific functional names.

---

## 1. Historical Stage Crosswalk

| Historical Stage ID | Canonical Scientific Module Name | Primary Implementation Directory | Locked Git Commit / Tag |
| :--- | :--- | :--- | :--- |
| **Step 1** | Project Initialization & Structure | Root, `config/`, `docs/architecture/` | `4cb87c0` |
| **Step 2** | System Architecture & Governance | `docs/architecture/SYSTEM_ARCHITECTURE.md` | `644e3f3` |
| **Step 3** | Metadata Contract & Schema | `config/metadata_schema.yaml`, `docs/data_dictionary/` | `6fb6dd8` |
| **Step 4** | Schema-Driven Metadata Validation | `src/metadata/`, `tests/metadata/` | `455265c` |
| **Step 5** | Ingestion & Modality-Aware Routing | `src/ingestion/`, `src/routing/` | `1bf8ddb` |
| **Step 6A** | Bulk Raw-Count Quality Control | `src/qc/`, `tests/qc/` | `3a8a7ef` |
| **Step 6B** | Canonical Gene Harmonization v1 | `src/gene_harmonization/`, `tests/gene_harmonization/` | `b904a41` |
| **Step 6C** | Bulk RNA-seq Normalization v1 | `src/normalization/`, `tests/normalization/` | `5ce67cc` |
| **Step 7A** | Treatment-versus-Matched-Control Contrast Builder v1 | `src/response_builder/`, `tests/response_builder/` | `cd954f4` |
| **Step 7B** | Bulk Differential Expression v1 | `src/differential_expression/`, `tests/differential_expression/` | `f1a7ad8` |
| **Web UI** | DNT Organoid AI Web Interface & Hardening | `app.py`, `src/app/`, `tests/app/` | `ca26d34` |
| **Stage 8C** | **Pathway and Signaling Features v1** | `src/features/pathway/`, `tests/features/pathway/` | **`stage8c-v1-frozen`** (`7a69c25`) |
| **Handoff v1** | Repository Organization & Handoff Documentation | `docs/project/`, root `README.md` | Current |

---

## 2. Chronological Milestone Ledger

### Milestone 1: Data Intake Foundation (Steps 1–5)
* Initialized DNT Organoid AI repository structure, data tiers, and architectural governance.
* Defined the canonical human brain-organoid metadata contract ([`config/metadata_schema.yaml`](../../config/metadata_schema.yaml)).
* Implemented schema-driven validation and asset routing for raw counts, FASTQs, UMIs, and single-cell objects.

### Milestone 2: Technical Standardization & Normalization (Steps 6A–6C)
* **Step 6A**: Implemented raw count QC with library size, zero-count fractions, and Spearman replicate correlation.
* **Step 6B**: Built automated Ensembl 112 + HGNC reference generator (`scripts/build_gene_reference.py`) and harmonizer preserving unmapped genes.
* **Step 6C**: Built official DESeq2 size-factor estimation bridge (`scripts/r/estimate_size_factors.R`) across metadata-defined cohorts.

### Milestone 3: Molecular Contrast & Differential Expression (Steps 7A–7B)
* **Step 7A**: Implemented `build_response_contrasts()` strictly pairing treatment organoids with vehicle controls sharing developmental age, cell line, and protocol ($\ge 2$ biological replicates per arm).
* **Step 7B**: Built live DESeq2 Wald test R bridge (`scripts/r/run_deseq2_contrast.R`) with gene-universe conservation auditing across all 60,207 canonical genes. Merged to `main` at `f1a7ad8`.

### Milestone 4: Streamlit Web Interface Hardening (Web App)
* Implemented Streamlit web interface (`app.py`, `src/app/`) supporting single-matrix uploads, manual sample assignments, official GEO metadata resolution (GSE232218), and interactive Volcano/MA visualization. Hardened at `ca26d34`.

### Milestone 5: Pathway & Signaling Feature Foundation (Stage 8C — Frozen)
* Implemented preranked GSEA for Gene Ontology (GO BP) and KEGG pathways using clusterProfiler 4.14.6.
* Established exact statistical universe preservation: 32,140 statistically eligible genes; 7,122 non-Entrez genes retained (`future_feature_eligibility=True`); 2 duplicate Entrez mappings resolved deterministically via `MAX_ABS_WALD`.
* Evaluated 4,893 GO BP pathways (1,151 significant) and 348 KEGG pathways (129 significant).
* Formally frozen with immutable annotated release tag: [`stage8c-v1-frozen`](https://github.com/bioprintingLabsInc/DNT-Organoid-AI/releases/tag/stage8c-v1-frozen) at commit `7a69c259b9be1ea9cff55dd76b8db365df105b9f`.

### Milestone 6: Repository Organization & Handoff Documentation (Current)
* Created comprehensive `docs/project/` documentation suite and standardized scientific module nomenclature.
* Mapped failure domains, runbooks, and lineages for new researchers and autonomous agents with zero runtime code modifications.
