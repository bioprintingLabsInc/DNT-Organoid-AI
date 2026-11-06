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
   Every gene evaluated in upstream differential expression (Stage 8A.1 / Stage 8A.2) is guaranteed representation in the evidence dataset. Genes not belonging to evaluated pathways or leading edges are preserved with explicit qualification (`GenePathwayMembershipType.CONFIRMED_NON_MEMBER` or `UNKNOWN_OR_UNAVAILABLE`, `PathwayEvidenceStatus.NOT_ANNOTATED` or `PathwayEvidenceStatus.NOT_AVAILABLE`). No arbitrary $p$-value or magnitude cutoffs are applied.
2. **Four-State Biological Gene–Pathway Membership:**
   The evidence layer strictly distinguishes four mutually exclusive membership states:
   - `LEADING_EDGE_MEMBER`: Verified core driver in preranked GSEA.
   - `NON_LEADING_EDGE_MEMBER`: Verified member of the complete pathway gene set outside the leading edge.
   - `CONFIRMED_NON_MEMBER`: A complete versioned reference confirms the gene is NOT a member of the pathway.
   - `UNKNOWN_OR_UNAVAILABLE`: A complete versioned gene-set reference was unavailable or incomplete, or pathway analysis was not performed.
3. **Decoupling of Pathway Significance from Membership:**
   Pathway significance (`is_pathway_significant: bool`, based on FDR < 0.05) is kept completely separate from topological gene membership (`gene_membership_type`).
4. **Decoupling of Perturbation from Direct Binding / Causality:**
   A transcriptomic fold change represents an *observed cellular perturbation*, not direct physical binding between a chemical ligand and a target protein. All records enforce `is_direct_target_interaction = False` and `causality_claim = "OBSERVED_TRANSCRIPTOMIC_RESPONSE"` unless verified direct assay evidence is provided.
5. **Decoupling of Gene Sign from Pathway Direction:**
   An individual gene's fold change sign ($\log_2\text{FC} > 0$ vs $< 0$) is mathematically decoupled from the pathway's enrichment direction (GSEA positive vs negative NES).
6. **Strict Contrast-Specific Integrity:**
   Evidence records strictly link genes and pathways to the exact study, chemical, treatment condition, control condition, and contrast. Cross-contrast mismatch between gene and pathway identities triggers immediate rejection rather than silent cross-joining.
7. **Universal, Metadata-Driven Architecture:**
   The module contains zero hardcoded chemical names, doses, durations, ages, or gene lists. Missing experimental metadata remains explicitly empty or unmeasured; the pipeline never invents sample identifiers, biological replicate IDs, or assumed replicate counts.
8. **Strict Leakage & Circularity Prevention:**
   The evidence layer contains zero references to regulatory DNT outcome labels (`src/labels/` is never imported).

---

## 2. Package Architecture

The Stage 8D implementation is located in `src/features/compound_pathway/`:

| Module | Purpose |
| :--- | :--- |
| [`models.py`](file:///src/features/compound_pathway/models.py) | Strongly-typed dataclasses: `ChemicalExposureContext`, `ChemicalGeneRecord`, `ChemicalPathwayRecord`, `ChemicalGenePathwayEvidenceRecord`, `ChemicalGenePathwayEvidenceDataset`, enums (`EvidenceType`, `PathwayEvidenceStatus`, `GenePathwayMembershipType`). |
| [`builder.py`](file:///src/features/compound_pathway/builder.py) | Lossless context extraction, deterministic evidence ID derivation (`_derive_evidence_id`), four-state membership evaluation (`determine_gene_pathway_membership`), contrast-specific linking, and high-performance inverted indexing. |
| [`serialization.py`](file:///src/features/compound_pathway/serialization.py) | Governed 43-column tabular export (`CANONICAL_EVIDENCE_COLUMNS`), DataFrame conversion, and deterministic CSV serialization. |
| [`qc.py`](file:///src/features/compound_pathway/qc.py) | Invariant auditing per study, contrast, and gene (`audit_evidence_dataset_qc`): verifies gene universe conservation, exact numerical parity against Stage 8A, metadata fidelity, duplicate evidence ID detection, cross-contrast integrity, and zero DNT label leakage. QC failure strictly blocks validation. |
| [`runner.py`](file:///src/features/compound_pathway/runner.py) | Generic execution entry point (`build_stage8d_evidence_pipeline`), lossless registry loader (`load_registry_from_csv`), robust boolean parser (`parse_bool_safely`), and pathway loader (`load_pathway_results_from_csv`). |

---

## 3. Governed Canonical Evidence Schema (43 Columns)

The output format is exported to a deterministic CSV with the following columns in exact sequence:

| # | Column Name | Type | Description |
| :- | :--- | :--- | :--- |
| 1 | `evidence_id` | str | Deterministic SHA-256 derived identifier (`ev_<hex16>`) |
| 2 | `study_id` | str | Governing study accession (e.g., `GSE232218`, `GSE187006`) |
| 3 | `contrast_id` | str | Treatment vs matched control contrast identifier |
| 4 | `independent_study_id` | str | Standardized independent study / study-family identifier |
| 5 | `compound_id` | str | Standard chemical identifier (e.g., `CID:3121`, `CID:2730`) |
| 6 | `compound_name` | str | Full compound name (e.g., `Valproic Acid`, `Chlorpyrifos`) |
| 7 | `treatment_condition_id` | str | Biological identifier of the treated experimental condition |
| 8 | `matched_control_condition_ids`| str | Semicolon-delimited vehicle control condition identifiers |
| 9 | `dose_value` | float | Tested chemical concentration |
| 10 | `dose_unit` | str | Concentration unit (e.g., `uM`, `mM`) |
| 11 | `exposure_duration` | float | Exposure time in specified units |
| 12 | `exposure_duration_unit` | str | Time unit (e.g., `h`) |
| 13 | `developmental_stage` | str | Organoid differentiation state at exposure start |
| 14 | `model_type` | str | Organoid model description |
| 15 | `canonical_ensembl_gene_id` | str | Primary Ensembl gene identifier (`ENSG...`) |
| 16 | `gene_symbol` | str | HGNC approved gene symbol |
| 17 | `base_mean` | float | Baseline expression level (DESeq2 baseMean or limma-trend mean) |
| 18 | `log2_fold_change` | float | Differential expression magnitude ($\log_2\text{FC}$) |
| 19 | `log2_fold_change_se` | float | Standard error of the $\log_2\text{FC}$ estimate |
| 20 | `wald_statistic` | float | Test statistic (Wald or moderated $t$) |
| 21 | `p_value` | float | Nominal differential expression significance |
| 22 | `adjusted_p_value` | float | Benjamini-Hochberg FDR |
| 23 | `gene_response_direction` | str | `UP`, `DOWN`, or `NO_CHANGE` |
| 24 | `gene_estimate_status` | str | `OK`, `LOW_COUNT`, `COOKS_OUTLIER`, etc. |
| 25 | `pathway_evidence_status` | str | `AVAILABLE`, `NOT_AVAILABLE`, `NOT_ANNOTATED`, `EVALUATION_INELIGIBLE` |
| 26 | `pathway_id` | str | Canonical pathway identifier (e.g. `GO:0007399`, `hsa04010`, or `None`) |
| 27 | `pathway_name` | str | Full pathway name |
| 28 | `source_database` | str | `GO_BP`, `KEGG`, or `None` |
| 29 | `database_version` | str | Version of underlying reference database |
| 30 | `enrichment_score` | float | GSEA raw Enrichment Score (ES) |
| 31 | `pathway_normalized_enrichment_score` | float | GSEA Normalized Enrichment Score (NES) |
| 32 | `pathway_p_value` | float | Nominal GSEA significance |
| 33 | `pathway_adjusted_p_value` | float | GSEA FDR |
| 34 | `pathway_size` | int | Total genes annotated in pathway gene set |
| 35 | `pathway_direction` | str | `POSITIVE`, `NEGATIVE`, or `None` |
| 36 | `is_pathway_significant` | bool | Boolean significance flag (e.g. FDR < 0.05) |
| 37 | `gene_membership_type` | str | `LEADING_EDGE_MEMBER`, `NON_LEADING_EDGE_MEMBER`, `CONFIRMED_NON_MEMBER`, `UNKNOWN_OR_UNAVAILABLE` |
| 38 | `is_leading_edge` | bool | Whether this specific gene is part of the GSEA leading edge |
| 39 | `evidence_type` | str | `OBSERVED_TRANSCRIPTOMIC_PERTURBATION`, `INACTIVE_OR_UNANNOTATED_PERTURBATION` |
| 40 | `is_direct_target_interaction`| bool | Strictly `False` for transcriptomic perturbation evidence |
| 41 | `causality_claim` | str | Strictly `"OBSERVED_TRANSCRIPTOMIC_RESPONSE"` |
| 42 | `scientific_limitations` | str | Explicit disclaimer of direct binding / directional pathway causality |
| 43 | `provenance_notes` | str | Analysis provenance, database sources, and version notes |

---

## 4. Real Data Validation & Reconciled Exposure Metadata

Stage 8D was executed across the validated human brain organoid cohorts via [`scripts/run_stage8d_evidence.py`](file:///scripts/run_stage8d_evidence.py).

> [!NOTE]
> **GSE232218 Metadata Reconciliation**:
> Exposure duration for GSE232218 is **48.0 hours** (`48h`), exactly matching the frozen source provenance ([`data/processed/GSE232218_VPA_200uM_provenance.json`](file:///data/processed/GSE232218_VPA_200uM_provenance.json)).

### Real Cohort Summary Table:

| Study & Regimen | Baseline Module | Input Genes | Conserved Genes | Total Evidence Records | Pathway Records | Unannotated Records | Leakage Violations | Numerical Mismatches | SHA-256 Checksum |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GSE232218**<br>(VPA 200 µM, **48h**) | Stage 8A.1 + Stage 8C (GO BP + KEGG) | 60,207 | 60,207 (100%) | 176,555 | 123,911 | 52,644 | 0 | 0 | `7cd636fe8653e81ac6d5850faaf77273b193cbddd6f6937bfe7a74764f849421` |
| **GSE187006 U1M**<br>(VPA 1 mM, 72h) | Stage 8A.2 limma-trend | 15,561 | 15,561 (100%) | 15,561 | 0 (`NOT_AVAILABLE`) | 15,561 | 0 | 0 | `5f741b3a27dcbce40d34d173c7105a7337e8d85467c8f5cb97b6f9e11c376f42` |
| **GSE187006 U2F**<br>(VPA 1 mM, 72h) | Stage 8A.2 limma-trend | 15,561 | 15,561 (100%) | 15,561 | 0 (`NOT_AVAILABLE`) | 15,561 | 0 | 0 | `dd4e65862f9aa56d60a671a2c4d15dec530065aa4ce0983045d202fc3ced980c` |

Manifest persisted at: [`data/processed/stage8d_artifacts_sha256.json`](file:///data/processed/stage8d_artifacts_sha256.json).

---

## 5. Verification & Test Execution

Run the complete Stage 8D test suite (15 tests):
```powershell
python -m pytest tests/features/compound_pathway/ -v -o pythonpath=". tests"
```
Output: **15/15 passed in 22.90s**.

Run the full feature regression suite (68 tests):
```powershell
python -m pytest tests/features/ -v -o pythonpath=". tests"
```
Output: **68/68 passed in 220.61s**.
