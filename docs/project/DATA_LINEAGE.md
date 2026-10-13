# Data Lineage and Provenance Map

This document details the transformations, schemas, storage locations, and provenance tracking across the DNT Organoid AI platform.

---

## 1. End-to-End Data Lifecycle

```text
Raw Data Sources (GEO/In-House)
       │
       ▼ [GSE232218_gene_expected_count.txt.gz]
1. Ingestion & Format Detection
       │  Artifact: RegisteredAsset (Format: raw_counts, Modality: bulk_rna_seq)
       ▼
2. Metadata Validation
       │  Artifact: ValidationResult (Schema: config/metadata_schema.yaml)
       ▼
3. Bulk Raw-Count QC
       │  Artifact: DatasetResult (Library sizes, zero fractions, Spearman correlations)
       ▼
4. Gene Harmonization
       │  Reference: data/references/ensembl_human_genes_v112.tsv
       │  Artifact: HarmonizedDataset (60,207 canonical Ensembl IDs mapped)
       ▼
5. Bulk RNA-seq Normalization
       │  Engine: DESeq2 ratio method (scripts/r/estimate_size_factors.R)
       │  Artifact: NormalizedDataset (Sample size factors + normalized counts)
       ▼
6. Treatment-Control Contrast Construction
       │  Criteria: Same study, experiment, organoid context, developmental age
       │  Artifact: MolecularResponseContrast (4 treated vs. 4 matched controls)
       ▼
7. Differential Molecular Response
       │  Engine: DESeq2 Wald test (scripts/r/run_deseq2_contrast.R)
       │  Artifact: DifferentialExpressionContrastResult (Wald, log2FC, FDR)
       ▼
8. Pathway & Signaling Features (Stage 8C)
       │  Engine: clusterProfiler 4.14.6 preranked GSEA (scripts/r/run_preranked_gsea.R)
       │  Artifacts:
       │    ├── data/processed/GSE232218_VPA_200uM_statistically_eligible_universe.csv
       │    ├── data/processed/GSE232218_VPA_200uM_retained_non_entrez_genes.csv
       │    ├── data/processed/GSE232218_VPA_200uM_GO_BP_results.csv
       │    ├── data/processed/GSE232218_VPA_200uM_KEGG_results.csv
       │    ├── data/processed/GSE232218_VPA_200uM_qc.json
       │    └── data/processed/GSE232218_VPA_200uM_provenance.json
```

---

## 2. Governed Reference Data Assets

All external knowledge bases and gene annotations are pinned and version-controlled locally:

1. **Canonical Human Gene Catalog**:
   * File: `data/references/ensembl_human_genes_v112.tsv`
   * Checksum (SHA-256): `3054aa734c15bb0a9cc64b67a8ead2a023e7dfd2cc0d0d6b7706372f4b13940e`
   * Release: Ensembl Release 112 (May 2024), Genome Build GRCh38.p14
   * Nomenclature Source: HGNC snapshot (2026-09-24)
   * Provenance: `data/references/ensembl_human_genes_v112_provenance.json`
   * Builder: `scripts/build_gene_reference.py`

2. **KEGG Pathway Database Offline Snapshot**:
   * File: `data/references/kegg_hsa_pathways.rds` & `kegg_hsa_pathways.json`
   * Checksum RDS (SHA-256): `A54298B77CE24D0059FB58CE0B3E22EB6837589D02FA4B2437E7545A02967C20`
   * Organism: *Homo sapiens* (`hsa`)
   * Pathway Count: 372 pathways (39,671 gene-pathway associations)
   * Provenance: `data/references/kegg_hsa_pathways_provenance.json`

3. **Gene Ontology Database (GO BP)**:
   * Engine: Bioconductor `org.Hs.eg.db` (v3.20.0) / `GO.db` (v3.20.0)
   * Key Type: `ENTREZID` (Bioconductor 3.20, R 4.4.3)
   * SQLite Checksum: `7b88d1cf970e8acef5da60b0535bbc1b53768381940d2b4fd538524f9b612d18`

---

## 3. Storage Layout & Data Tiers

* `data/raw/`: Immutable raw inputs (FASTQ files, vendor matrices, raw count tables).
* `data/interim/`: Intermediate serialized dataclasses during pipeline execution.
* `data/processed/`: Tabular feature tables, differential expression tables, and GSEA enrichment results.
* `data/model_ready/`: Frozen feature matrices and labels partitioned for training/validation.
* `data/manifests/`: YAML manifests tracking source asset checksums and acquisition dates.
* `outputs/`: Human-readable summary PDFs, Volcano plots, MA plots, and QC run reports.
