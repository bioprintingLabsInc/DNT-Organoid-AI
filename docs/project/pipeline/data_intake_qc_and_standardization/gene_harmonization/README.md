# Gene Harmonization

**Historical Identifier**: Step 6B  
**Validation Status**: **VALIDATED** (Commit: `b904a41`)  
**Umbrella Position**: Data Intake — Failure Domain 4

---

## 1. Scientific & Computational Purpose
Maps heterogeneous input gene identifiers (symbols, previous symbols, aliases, Entrez IDs, versioned Ensembl IDs) to canonical **Ensembl Release 112** IDs (`GRCh38.p14`). Preserves original identifiers, expression values, and unmapped genes without synthetic aggregation or silent deletion.

---

## 2. Inputs & Outputs
* **Inputs**: Count matrix (from QC step) + validated metadata + governed reference catalog.
* **Outputs**: `HarmonizedDataset` containing harmonized gene records (`HarmonizedGene`), mapping status (`CANONICAL`, `MAPPED_FROM_SYMBOL`, `UNMAPPED`), and collision audit trails.

---

## 3. Code, Config & Reference Map
* **Primary Code**:
  * [`src/gene_harmonization/type_detector.py`](../../../../../src/gene_harmonization/type_detector.py): Identifier type resolver.
  * [`src/gene_harmonization/reference.py`](../../../../../src/gene_harmonization/reference.py): In-memory Ensembl 112 lookup index.
  * [`src/gene_harmonization/harmonizer.py`](../../../../../src/gene_harmonization/harmonizer.py): Orchestration engine (`harmonize_bulk_counts`).
* **Configuration & Reference**:
  * [`config/gene_harmonization.yaml`](../../../../../config/gene_harmonization.yaml): Supported types and regex patterns.
  * [`data/references/ensembl_human_genes_v112.tsv`](../../../../../data/references/ensembl_human_genes_v112.tsv): Governed catalog.
  * [`scripts/build_gene_reference.py`](../../../../../scripts/build_gene_reference.py): Reference builder script.
* **Test Suite**:
  * [`tests/gene_harmonization/test_harmonizer.py`](../../../../../tests/gene_harmonization/test_harmonizer.py)
  * [`tests/gene_harmonization/test_reference.py`](../../../../../tests/gene_harmonization/test_reference.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: High unmapped gene rate or collision error**:
  * *First file to inspect*: [`src/gene_harmonization/type_detector.py`](../../../../../src/gene_harmonization/type_detector.py).
  * *Second file to inspect*: [`src/gene_harmonization/harmonizer.py`](../../../../../src/gene_harmonization/harmonizer.py).
  * *Common cause*: Input matrix contains non-human identifiers or unstripped transcript version tags.
