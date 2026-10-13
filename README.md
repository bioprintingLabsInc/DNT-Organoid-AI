# DNT Organoid AI

> **A research platform for developing and validating AI models that assess developmental-neurotoxicity-associated molecular responses in human brain organoids.**

👉 **New researchers, developers, and AI agents: Start with the [Quick Orientation Guide](docs/project/START_HERE.md).**

---

## 1. Project Purpose & Scope

DNT Organoid AI translates genome-wide transcriptomic profiling (bulk and single-cell RNA-seq) of human brain organoids exposed to chemical compounds into interpretable, leakage-safe toxicological features and predictive machine learning models.

* **Experimental Paradigm**: Organoid treatments are strictly compared against matched vehicle controls sharing identical developmental age and culture protocols ($\ge 2$ biological replicates per arm).
* **Zero Leakage**: Feature extraction operates strictly on the cellular molecular perturbation without using DNT outcome labels or expected biological answers.
* **Dual-Stack Execution**: High-performance Python 3.12 data orchestration coupled with official Bioconductor R 4.4.3 statistical engines (`DESeq2 1.46.0`, `clusterProfiler 4.14.6`).

---

## 2. Current Validated State & Frozen Baseline

* **Production Trunk (`main`)**: Data Intake, QC & Standardization, Treatment-Control Response, and Differential Molecular Response validated.
* **Web Interface**: Interactive Streamlit web interface (`app.py`) validated for manual sample assignment and GEO accession lookup (`ca26d34`).
* **Frozen Scientific Baseline**:
  * **Tag**: [`stage8c-v1-frozen`](https://github.com/bioprintingLabsInc/DNT-Organoid-AI/tree/stage8c-v1-frozen)
  * **Commit**: `7a69c259b9be1ea9cff55dd76b8db365df105b9f`
  * **Milestone**: **Pathway and Signaling Features v1** (GO BP & KEGG preranked GSEA).
  * **Accounting on GSE232218**: 60,207 evaluated genes; 32,140 statistically eligible genes; 7,122 non-Entrez valid genes retained; 25,016 final GSEA ranked genes; 4,893 GO BP pathways and 348 KEGG pathways tested.

---

## 3. Scientific Pipeline Hierarchy

```text
Data Intake, QC & Standardization
├── Ingestion & Format Detection
├── Metadata Validation
├── Quality Control & Preprocessing
├── Gene Harmonization (Ensembl 112)
└── Normalization (DESeq2 Size Factors)

Treatment-Control Molecular Response
├── Treatment-Control Response Construction
└── Differential Molecular Response (DESeq2 Wald Test)

Molecular Feature Engineering
├── Pathway & Signaling Features [FROZEN v1]
├── Cross-Study Gene Response Features [PLANNED]
├── Gene Coresponse & Network Features [PLANNED]
└── Compound-Gene-Pathway Relationships [PLANNED]

Evidence & Dataset Partitioning
├── Independent DNT Label Curation [PLANNED]
└── Compound/Study-Safe Dataset Splitting [PLANNED]

Predictive Modeling & Evaluation
├── DNT Signature Discovery (Training Partition Only) [PLANNED]
├── Model Training [PLANNED]
├── Probabilistic Calibration [PLANNED]
└── Untouched Internal Test Evaluation [PLANNED]
```
For the complete visual flowchart, data transition contracts, and leakage-safe rules, see the **[Pipeline Map](docs/project/PIPELINE_MAP.md)**.

---

## 4. Quickstart: Running Tests & Reproducing Results

### Run Automated Tests:
```bash
py -3.12 -m pytest -q
```
*Current test health: 327 passed, 13 subtests passed.*

### Run the Interactive Web Interface:
```bash
py -3.12 -m streamlit run app.py
```
*Web interface accessible at `http://localhost:8501`.*

### Reproduce Public Benchmark (GSE232218 Valproic Acid):
```bash
py -3.12 scripts/run_vpa_pilot_stage8c.py
```
> [!NOTE]
> **Environment & Portability Notice**: `scripts/run_vpa_pilot_stage8c.py` currently targets a local path on the development machine (`C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt.gz`). It is fully validated in the development environment. For running on other machines, see the **[Reproduction Runbook](docs/project/REPRODUCTION_RUNBOOK.md)** for exact file placement requirements until CLI path parameterization is implemented in a future code update.

---

## 5. Troubleshooting & Diagnostics

If an input dataset, test, or bridge fails, refer to the **[Practical Troubleshooting Guide](docs/project/TROUBLESHOOTING.md)** to diagnose across the 5 intake failure domains:
1. *File cannot be read / format mismatch* → [`src/routing/detector.py`](src/routing/detector.py)
2. *Metadata schema invalid* → [`src/metadata/validator.py`](src/metadata/validator.py)
3. *Raw count QC failed* → [`src/qc/validator.py`](src/qc/validator.py)
4. *Gene ID unmapped / collision* → [`src/gene_harmonization/harmonizer.py`](src/gene_harmonization/harmonizer.py)
5. *Normalization / size factor failed* → [`src/normalization/normalizer.py`](src/normalization/normalizer.py)

---

## 6. Regulatory Disclaimer

This repository is an academic and industrial research development platform. It is not a clinical, diagnostic, or regulatory safety decision system, and its outputs must not be used as such.
