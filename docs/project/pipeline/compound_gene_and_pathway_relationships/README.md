# Stage 8D: Chemical–Gene–Pathway Evidence Relationships

**Canonical Scientific Name:** Compound-Gene and Compound-Pathway Relationships  
**Historical Identifier:** Stage 8D  
**Status:** **VALIDATED & IMPLEMENTED**  
**Package:** `src/features/compound_pathway/`  
**Test Suite:** `tests/features/compound_pathway/test_evidence_relationships.py`  
**Orchestration Script:** `scripts/run_stage8d_evidence.py`

---

## 1. Scientific Objective & Governing Design

Stage 8D synthesizes structured, traceable, and biologically faithful evidence relationships connecting experimental chemical exposures to observed molecular perturbations in human brain organoids:

$$\text{Chemical Exposure} \longrightarrow \text{Observed Gene Response} \longrightarrow \text{Associated Biological Pathway}$$

### Core Scientific Principles:
1. **Law of Conservation of the Gene Universe:**
   Every gene evaluated in upstream differential expression (Stage 8A.1) is guaranteed representation in the evidence dataset. Genes not belonging to enriched pathways or leading edges are preserved with explicit qualification (`GenePathwayMembershipType.NOT_MEMBER`, `PathwayEvidenceStatus.NOT_ANNOTATED` or `PathwayEvidenceStatus.NOT_AVAILABLE`). No arbitrary $p$-value or magnitude cutoffs are applied.
2. **Decoupling of Observed Perturbation from Direct Binding / Causality:**
   A transcriptomic fold change indicates an *observed biological response*, not direct physical chemical–target binding. All records enforce `is_direct_target_interaction = False` and `causality_claim = "OBSERVED_TRANSCRIPTOMIC_RESPONSE"` unless direct biophysical assay proof is explicitly supplied.
3. **Decoupling of Gene Sign from Pathway Direction:**
   The sign of an individual gene's fold change ($\log_2\text{FC} > 0$ vs $< 0$) is mathematically decoupled from the pathway's enrichment direction (GSEA positive vs negative NES). A downregulated gene can participate as a leading-edge driver in an activated pathway, and vice-versa.
4. **Pathway Availability Gating:**
   Contrasts lacking valid empirical GSEA pathway runs (e.g. GSE187006 U1M and U2F) gracefully degrade to `PathwayEvidenceStatus.NOT_AVAILABLE` with explicit scientific limitation notes, preserving 100% of the gene universe without fabricating associations.
5. **Universal, Metadata-Driven Architecture:**
   The module contains zero hardcoded chemical names, doses, durations, ages, or gene lists. It generalizes across any chemical entity (e.g. VPA, Chlorpyrifos, Rotenone) and organoid regimen.
6. **Strict Leakage & Circularity Prevention:**
   The evidence layer contains zero references to regulatory DNT outcome labels (`src/labels/` is never imported). External toxicological mechanisms are isolated from input features.

---

## 2. Package Architecture

The Stage 8D implementation is located in `src/features/compound_pathway/`:

| Module | Purpose |
| :--- | :--- |
| [`models.py`](file:///src/features/compound_pathway/models.py) | Strongly-typed dataclasses: `ChemicalExposureContext`, `ChemicalGeneRecord`, `ChemicalPathwayRecord`, `ChemicalGenePathwayEvidenceRecord`, `ChemicalGenePathwayEvidenceDataset`, enums (`EvidenceType`, `PathwayEvidenceStatus`, `GenePathwayMembershipType`). |
| [`builder.py`](file:///src/features/compound_pathway/builder.py) | Context extraction, deterministic evidence ID derivation (`_derive_evidence_id`), chemical-gene construction, pathway linking, and gene universe reconciliation. |
| [`serialization.py`](file:///src/features/compound_pathway/serialization.py) | Governed 37-column tabular export (`CANONICAL_EVIDENCE_COLUMNS`), DataFrame conversion, and deterministic CSV serialization. |
| [`qc.py`](file:///src/features/compound_pathway/qc.py) | Invariant auditing (`audit_evidence_dataset_qc`): verifies gene universe conservation, exact numerical identity against Stage 8A.1, zero DNT label leakage, and absence of speculative causality claims. |
| [`runner.py`](file:///src/features/compound_pathway/runner.py) | Dataset-independent generic pipeline (`build_stage8d_evidence_pipeline`), registry loader (`load_registry_from_csv`), and pathway loader (`load_pathway_results_from_csv`). |

---

## 3. Canonical Evidence Schema (37 Columns)

The output format is exported to a deterministic CSV with the following columns:

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `evidence_id` | str | Deterministic SHA-256 derived identifier for the evidence link |
| `study_id` | str | GEO or BioProject accession (e.g., `GSE232218`, `GSE187006`) |
| `contrast_id` | str | Treatment vs matched control contrast identifier |
| `compound_id` | str | Standard chemical identifier (e.g., `CID:3121`, `CID:2730`) |
| `compound_name` | str | Full compound name (e.g., `Valproic Acid`, `Chlorpyrifos`) |
| `dose_value` | float | Tested chemical concentration |
| `dose_unit` | str | Concentration unit (e.g., `uM`, `mM`) |
| `dose_original` | str | Original reported dosage string |
| `exposure_duration` | float | Exposure time in specified units |
| `exposure_duration_unit` | str | Time unit (e.g., `h`) |
| `developmental_stage` | str | Organoid differentiation state at exposure start |
| `model_type` | str | Organoid model description |
| `canonical_ensembl_gene_id` | str | Primary Ensembl gene identifier (`ENSG...`) |
| `gene_symbol` | str | HGNC approved gene symbol |
| `base_mean` | float | Baseline expression level (DESeq2 baseMean or limma-trend mean) |
| `log2_fold_change` | float | Differential expression magnitude ($\log_2\text{FC}$) |
| `log2_fold_change_se` | float | Standard error of the $\log_2\text{FC}$ estimate |
| `wald_statistic` | float | Test statistic (Wald or moderated $t$) |
| `p_value` | float | Nominal differential expression significance |
| `adjusted_p_value` | float | Benjamini-Hochberg FDR |
| `gene_response_direction` | str | `UP`, `DOWN`, or `NO_CHANGE` |
| `gene_estimate_status` | str | `OK`, `LOW_COUNT`, `OUTLIER`, etc. |
| `pathway_id` | str | Canonical pathway identifier (e.g. `GO:0007399`, `hsa04010`, or `NONE`) |
| `pathway_name` | str | Full pathway name |
| `source_database` | str | `GO_BP`, `KEGG`, or `NONE` |
| `enrichment_score` | float | GSEA enrichment score (ES) |
| `normalized_enrichment_score` | float | GSEA normalized enrichment score (NES) |
| `pathway_p_value` | float | Nominal GSEA significance |
| `pathway_adjusted_p_value` | float | GSEA FDR |
| `pathway_size` | int | Total genes annotated in pathway gene set |
| `pathway_direction` | str | `POSITIVE`, `NEGATIVE`, or `NONE` |
| `is_leading_edge` | bool | Whether this specific gene is part of the GSEA leading edge |
| `membership_type` | str | `LEADING_EDGE`, `PATHWAY_MEMBER`, or `NOT_MEMBER` |
| `evidence_type` | str | Evidence classification category |
| `pathway_evidence_status` | str | `AVAILABLE`, `NOT_AVAILABLE`, `NOT_ANNOTATED` |
| `is_direct_target_interaction`| bool | Strictly `False` for transcriptomic perturbation evidence |
| `causality_claim` | str | Strictly `"OBSERVED_TRANSCRIPTOMIC_RESPONSE"` |

---

## 4. Real Data Validation & Artifact Checksums

Stage 8D was executed across the validated human brain organoid cohorts via [`scripts/run_stage8d_evidence.py`](file:///scripts/run_stage8d_evidence.py).

### Real Cohort Summary Table:

| Study & Regimen | Baseline Module | Input Genes | Conserved Genes | Total Evidence Records | Pathway Records | Unannotated Records | Leakage Violations | Numerical Mismatches | SHA-256 Checksum |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GSE232218**<br>(VPA 200 µM, 8h) | Stage 8A.1 + Stage 8C (GO BP + KEGG) | 60,207 | 60,207 (100%) | 176,555 | 123,911 | 52,644 | 0 | 0 | `03f916ba39338f68d8b16e934d3c1aa992e1fb8839938103e24a3e99be81406f` |
| **GSE187006 U1M**<br>(VPA 1 mM, 72h) | Stage 8A.2 limma-trend | 15,561 | 15,561 (100%) | 15,561 | 0 (`NOT_AVAILABLE`) | 15,561 | 0 | 0 | `57be587cf66e2b32bf42ac007e4f0695d943aefd3514e2493096fe82bd319edb` |
| **GSE187006 U2F**<br>(VPA 1 mM, 72h) | Stage 8A.2 limma-trend | 15,561 | 15,561 (100%) | 15,561 | 0 (`NOT_AVAILABLE`) | 15,561 | 0 | 0 | `159d759854a8998d3f174e2b1068795d182b2a9027d3d53a78dca5d11709375f` |

Manifest persisted at: [`data/processed/stage8d_artifacts_sha256.json`](file:///data/processed/stage8d_artifacts_sha256.json).

---

## 5. Verification & Test Execution

Run the complete Stage 8D test suite:
```powershell
python -m pytest tests/features/compound_pathway/ -v -o pythonpath=". tests"
```
Output: **8/8 passed in 19.90s**.

Run the full feature regression suite:
```powershell
python -m pytest tests/features/ -v -o pythonpath=". tests"
```
Output: **61/61 passed in 133.40s**.
