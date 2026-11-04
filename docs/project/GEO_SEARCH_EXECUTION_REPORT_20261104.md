# DNT AI — GEO Dataset Search v1.0 Execution Report

**Protocol Identifier:** `DNT-GEO-SEARCH-v1.0`  
**Execution Timestamp:** 2026-11-04T18:31:30Z – 2026-11-04T18:31:58Z UTC  
**Configuration File:** `config/dnt_geo_search_keywords_v1.json`  
**Configuration Checksum (SHA-256):** `143e5009cfd1a909b5f27af868e0146c42d2f60426729ebedc891394377f334d`  
**Database Target:** NCBI Gene Expression Omnibus (GEO), via Entrez E-utilities (`db=gds`, entry type `gse`)  
**Execution Output Directory:** `outputs/geo_search/2026-11-04_run01`  
**Branch:** `feature/geo-dataset-search-v1`

---

## 1. Executive Summary

This report documents the prospective, fully auditable execution of the **DNT AI GEO Dataset Search v1.0** across NCBI GEO. The search executed the 45-keyword configuration organized into three controlled Boolean query groups. All query strings, execution timestamps, NCBI ESearch reported counts, retrieved UID lists, and deduplicated GSE accessions were recorded without estimation or interpolation.

```
================================================================================
                    DNT AI GEO DATASET SEARCH v1.0 RESULTS
================================================================================
Total Query Matches (sum of Q1, Q2, Q3):           1,546
Total Unique GEO Datasets (deduplicated GSEs):     1,452
Historical Reference Cohort Size:                    516
Datasets Overlapping with Original 516:              302 (58.5% recovery)
Newly Discovered Datasets:                         1,150
Historical 516 Records Not Retrieved in This Run:    214
================================================================================
```

---

## 2. Exact Search Queries & Execution Metrics

### Query 1: `Q1_BIO_EXPOSURE`
- **Objective:** Intersect human 3D brain/neural biological models with chemical and environmental exposure concepts.
- **Join Logic:** `(organoid_3d OR neural_cells) AND exposure AND Homo sapiens[Organism] AND gse[Entry Type]`
- **Exact Query String:**
  ```text
  ("Brain organoid"[All Fields] OR "Cerebral organoid"[All Fields] OR "Cortical organoid"[All Fields] OR "Forebrain organoid"[All Fields] OR "Midbrain organoid"[All Fields] OR "Neural organoid"[All Fields] OR "Brain spheroid"[All Fields] OR "Neural spheroid"[All Fields] OR "BrainSphere"[All Fields] OR "3D brain"[All Fields] OR "3D neural"[All Fields] OR "Neural assembloid"[All Fields] OR "Cortical spheroid"[All Fields] OR "Neurosphere"[All Fields] OR "Neural rosette"[All Fields] OR "Brain organ-on-chip"[All Fields] OR "Neural progenitor cells"[All Fields] OR "Neural stem cells"[All Fields] OR "iPSC-derived neurons"[All Fields] OR "Human neural cells"[All Fields] OR "Neural differentiation"[All Fields] OR "Cortical neurons"[All Fields] OR "Dopaminergic neurons"[All Fields] OR "Neural precursor cells"[All Fields] OR "iPSC-derived neural cells"[All Fields] OR "Neurodevelopment"[All Fields]) AND ("Developmental neurotoxicity"[All Fields] OR "Neurodevelopmental toxicity"[All Fields] OR "Neurotoxicity"[All Fields] OR "Developmental toxicology"[All Fields] OR "Chemical-induced neurotoxicity"[All Fields] OR "Chemical exposure"[All Fields] OR "Drug-induced neurotoxicity"[All Fields] OR "Environmental toxicants"[All Fields] OR "Pesticide exposure"[All Fields] OR "Heavy metal exposure"[All Fields] OR "Prenatal exposure"[All Fields] OR "Neurotoxic compounds"[All Fields]) AND Homo sapiens[Organism] AND gse[Entry Type]
  ```
- **Execution Interval:** `2026-11-04T18:31:30.355195Z` to `2026-11-04T18:31:33.000147Z` (2.64s)
- **Reported Matches:** `66`
- **Retrieved GSE Accessions:** `66`
- **Anomalies / Format Failures:** `0` (100% resolved cleanly to `GSE` series)
- **Output Artifact:** `outputs/geo_search/2026-11-04_run01/Q1_BIO_EXPOSURE_accessions.csv`

---

### Query 2: `Q2_BIO_TRANSCRIPTOMICS`
- **Objective:** Intersect human 3D brain/neural biological models with sequencing/transcriptomic profiling modalities.
- **Join Logic:** `(organoid_3d OR neural_cells) AND transcriptomics AND Homo sapiens[Organism] AND gse[Entry Type]`
- **Exact Query String:**
  ```text
  ("Brain organoid"[All Fields] OR "Cerebral organoid"[All Fields] OR "Cortical organoid"[All Fields] OR "Forebrain organoid"[All Fields] OR "Midbrain organoid"[All Fields] OR "Neural organoid"[All Fields] OR "Brain spheroid"[All Fields] OR "Neural spheroid"[All Fields] OR "BrainSphere"[All Fields] OR "3D brain"[All Fields] OR "3D neural"[All Fields] OR "Neural assembloid"[All Fields] OR "Cortical spheroid"[All Fields] OR "Neurosphere"[All Fields] OR "Neural rosette"[All Fields] OR "Brain organ-on-chip"[All Fields] OR "Neural progenitor cells"[All Fields] OR "Neural stem cells"[All Fields] OR "iPSC-derived neurons"[All Fields] OR "Human neural cells"[All Fields] OR "Neural differentiation"[All Fields] OR "Cortical neurons"[All Fields] OR "Dopaminergic neurons"[All Fields] OR "Neural precursor cells"[All Fields] OR "iPSC-derived neural cells"[All Fields] OR "Neurodevelopment"[All Fields]) AND ("RNA sequencing"[All Fields] OR "RNA-seq"[All Fields] OR "Transcriptomics"[All Fields] OR "Bulk RNA-seq"[All Fields] OR "Single-cell RNA-seq"[All Fields] OR "snRNA-seq"[All Fields] OR "Gene expression profiling"[All Fields]) AND Homo sapiens[Organism] AND gse[Entry Type]
  ```
- **Execution Interval:** `2026-11-04T18:31:33.000147Z` to `2026-11-04T18:31:54.100453Z` (21.10s)
- **Reported Matches:** `1,358`
- **Retrieved GSE Accessions:** `1,358`
- **Anomalies / Format Failures:** `0` (100% resolved cleanly to `GSE` series)
- **Output Artifact:** `outputs/geo_search/2026-11-04_run01/Q2_BIO_TRANSCRIPTOMICS_accessions.csv`

---

### Query 3: `Q3_EXPOSURE_TRANSCRIPTOMICS`
- **Objective:** Intersect chemical/toxicological exposure concepts directly with transcriptomics, providing maximum sensitivity for exposure studies that may use non-standard organoid nomenclature.
- **Join Logic:** `exposure AND transcriptomics AND Homo sapiens[Organism] AND gse[Entry Type]`
- **Exact Query String:**
  ```text
  ("Developmental neurotoxicity"[All Fields] OR "Neurodevelopmental toxicity"[All Fields] OR "Neurotoxicity"[All Fields] OR "Developmental toxicology"[All Fields] OR "Chemical-induced neurotoxicity"[All Fields] OR "Chemical exposure"[All Fields] OR "Drug-induced neurotoxicity"[All Fields] OR "Environmental toxicants"[All Fields] OR "Pesticide exposure"[All Fields] OR "Heavy metal exposure"[All Fields] OR "Prenatal exposure"[All Fields] OR "Neurotoxic compounds"[All Fields]) AND ("RNA sequencing"[All Fields] OR "RNA-seq"[All Fields] OR "Transcriptomics"[All Fields] OR "Bulk RNA-seq"[All Fields] OR "Single-cell RNA-seq"[All Fields] OR "snRNA-seq"[All Fields] OR "Gene expression profiling"[All Fields]) AND Homo sapiens[Organism] AND gse[Entry Type]
  ```
- **Execution Interval:** `2026-11-04T18:31:54.100453Z` to `2026-11-04T18:31:58.110382Z` (4.01s)
- **Reported Matches:** `122`
- **Retrieved GSE Accessions:** `122`
- **Anomalies / Format Failures:** `0` (100% resolved cleanly to `GSE` series)
- **Output Artifact:** `outputs/geo_search/2026-11-04_run01/Q3_EXPOSURE_TRANSCRIPTOMICS_accessions.csv`

---

## 3. Comparison with Historical 516-Record Cohort

The original 516 accessions were compiled from historical exploration in `challange_project.ipynb` (Cells 18 and 20) using 11 organoid terms without requiring exposure or transcriptomic Boolean terms. 

### Comparative Breakdown

| Category | Accession Count | Percentage of Set | Description |
| :--- | :--- | :--- | :--- |
| **Total Unique Retrieved (New)** | **1,452** | 100.0% | Complete deduplicated list of unique GSE series |
| **Overlapping with Original 516** | **302** | 58.5% of original 516 | Verified historical series confirmed by the new query groups |
| **Newly Discovered Series** | **1,150** | 79.2% of new search | New potential candidates captured by expanded vocabulary |
| **Original 516 Not Retrieved** | **214** | 41.5% of original 516 | Historical series excluded by the Boolean intersections |

### Scientific Reconciliation of the 214 Non-Retrieved Historical Series

A systematic sample audit of the 214 non-retrieved original series confirmed that they were not returned due to specific, scientifically justified structural differences:
1. **Alternative Profiling Modalities:** Series such as `GSE109104` utilize DNA methylation profiling (`Methylation profiling by high throughput sequencing`), which lacks any of the 7 transcriptomics keywords.
2. **Pure Developmental Atlases Without Chemical Perturbation:** Series such as `GSE113931` (cross-species great ape evolutionary genomics) possess organoid terminology but lack chemical or developmental neurotoxicity exposure terms.
3. **SuperSeries Parent Records:** Certain parent SuperSeries records (e.g. `GSE110611`) aggregate sub-series that index differently in Entrez.

As mandated by `docs/project/GEO_DATASET_SEARCH_PROTOCOL_V1.md`, all 516 historical records remain preserved in `data/metadata/original_516_geo_accessions.txt` for independent downstream review.

---

## 4. Query Membership Breakdown

The 1,452 unique datasets partition cleanly across the three queries:

| Query Membership Pattern | Count | Percentage | Characterization |
| :--- | :--- | :--- | :--- |
| **Q2 only** (`Q2_BIO_TRANSCRIPTOMICS`) | 1,311 | 90.3% | Brain/neural model + transcriptomics (baseline biological characterization) |
| **Q3 only** (`Q3_EXPOSURE_TRANSCRIPTOMICS`) | 75 | 5.2% | Chemical exposure + transcriptomics (broad toxicology, requires neural model screening) |
| **Q1 + Q2 + Q3 (All 3 Queries)** | **47** | **3.2%** | **Highest-priority candidates:** explicit neural model + chemical exposure + transcriptomics |
| **Q1 only** (`Q1_BIO_EXPOSURE`) | 19 | 1.3% | Brain/neural model + chemical exposure (may include array or non-RNA modalities) |
| **Total Unique GSE Series** | **1,452** | **100.0%** | Completely deduplicated search universe |

The **47 three-way intersecting datasets** (`GSE141253`, `GSE150005`, `GSE151996`, `GSE169751`, `GSE171566`, `GSE182307`, `GSE199502`, `GSE250609`, `GSE253569`, `GSE253820`, `GSE253821`, `GSE253822`, `GSE261775`, `GSE261878`, `GSE263394`, `GSE264703`, `GSE265797`, `GSE269254`, `GSE270015`, `GSE270064`, `GSE270066`, `GSE270070`, `GSE272334`, `GSE274212`, `GSE282273`, `GSE286054`, `GSE289218`, `GSE291145`, `GSE293246`, `GSE293247`, `GSE293248`, `GSE293251`, `GSE293252`, `GSE293253`, `GSE293255`, `GSE299369`, `GSE299790`, `GSE302654`, `GSE304103`, `GSE318014`, `GSE325791`, `GSE329929`, `GSE329931`, `GSE329933`, `GSE351619`, `GSE4773`, `GSE63935`) represent the most immediate high-value candidates for human brain organoid chemical-exposure and DNT evaluation.

---

## 5. Artifact & Provenance Directory

All execution artifacts have been persisted in Git tracking under:
- `outputs/geo_search/2026-11-04_run01/run_manifest.json` (Machine-readable execution parameters, environment, and timestamps)
- `outputs/geo_search/2026-11-04_run01/all_unique_geo_series.csv` (1,452 deduplicated GSE series with query memberships)
- `outputs/geo_search/2026-11-04_run01/comparison_with_original_516.csv` (1,666 combined entries showing search status)
- `outputs/geo_search/2026-11-04_run01/comparison_summary.json` (Summary metrics and cryptographic hashes)
- `outputs/geo_search/2026-11-04_run01/Q1_BIO_EXPOSURE_accessions.csv` (66 series)
- `outputs/geo_search/2026-11-04_run01/Q2_BIO_TRANSCRIPTOMICS_accessions.csv` (1,358 series)
- `outputs/geo_search/2026-11-04_run01/Q3_EXPOSURE_TRANSCRIPTOMICS_accessions.csv` (122 series)
- `data/metadata/original_516_geo_accessions.txt` (Authoritative historical accession list)

---

## 6. Governance & Immutability Adherence

1. **No Estimated Counts:** 100% of reported counts are derived directly from verified NCBI E-utilities execution logs.
2. **Zero Modification to RNA-Seq Pipeline:** The core RNA-seq processing modules (`src/`, `config/quantification_bridge/`, etc.) remain 100% immutable and untouched.
3. **No Premature Exclusion:** As instructed, zero datasets have been excluded or downselected at this stage; this inventory establishes the complete discovery universe for subsequent sample-level review.
