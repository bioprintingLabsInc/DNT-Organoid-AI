# Welcome to DNT Organoid AI: Quick Orientation

Welcome to the **DNT Organoid AI** codebase. This guide is designed to orient new scientists, software engineers, and AI agents within 5 minutes.

---

## 1. Executive Summary

* **Scientific Mission**: Developing and validating an interpretable, leakage-safe machine learning platform that evaluates developmental neurotoxicity (DNT) molecular responses in human brain organoids from bulk and single-cell RNA-sequencing data.
* **Core Paradigm**: The platform evaluates **treatment-versus-matched-control molecular perturbations** without assuming prior biological answers or leaking DNT outcome labels into feature extraction.
* **Current Validated Milestone**: Stages 1 through 7B (Data Intake, Normalization, Contrast Construction, Differential Expression) and Stage 8C (Pathway and Signaling Feature Foundation) are fully implemented and verified against public human cortical organoid data (**GSE232218**, Valproic Acid 200 µM).
* **Current Frozen Baseline**: Tag [`stage8c-v1-frozen`](https://github.com/bioprintingLabsInc/DNT-Organoid-AI/releases/tag/stage8c-v1-frozen) at commit `7a69c259b9be1ea9cff55dd76b8db365df105b9f`.

---

## 2. Where Things Live

```text
DNT-Organoid-AI/
├── app.py                     # Streamlit web interface entry point
├── config/                    # Governed YAML schemas and QC rules
├── data/                      # Manifests, metadata, processed universes, and references
│   └── references/            # Ensembl Release 112 catalog and offline KEGG/GO assets
├── docs/                      # Architectural specifications and project runbooks
│   └── project/               # Primary developer & scientist handoff system
│       └── pipeline/          # Scientific module documentation packages
├── scripts/                   # Reference builders and pilot execution runners
│   └── r/                     # Headless R engines (DESeq2, clusterProfiler, renv.lock)
├── src/                       # Production Python packages
│   ├── app/                   # Web interface backend orchestration and GEO resolver
│   ├── ingestion/             # Asset registration
│   ├── metadata/              # Schema validation
│   ├── routing/               # Format and modality detection
│   ├── qc/                    # Bulk raw-count quality control
│   ├── gene_harmonization/    # Canonical Ensembl ID harmonization
│   ├── normalization/         # DESeq2 size-factor normalization
│   ├── response_builder/      # Treatment-vs-matched-control contrasts
│   ├── differential_expression/ # DESeq2 Wald differential expression
│   └── features/pathway/      # GO BP and KEGG preranked GSEA [FROZEN v1]
└── tests/                     # 327 passing automated tests + 13 subtests
```

---

## 3. Five-Minute Quickstart

### Step 1: Environment Setup
Ensure Python 3.12 and R 4.4.3 (with Bioconductor 3.20) are installed:
```bash
# Install Python dependencies
py -3.12 -m pip install -r requirements.txt
```

### Step 2: Run the Automated Test Suite
Verify that all 327 unit and integration tests pass:
```bash
py -3.12 -m pytest -q
```
*Expected Output*: `327 passed, 13 subtests passed in ~6 minutes`.

### Step 3: Run the Streamlit Web Interface
Launch the interactive web interface locally:
```bash
py -3.12 -m streamlit run app.py
```
*Access in browser*: `http://localhost:8501`.

### Step 4: Reproduce the Pilot Dataset (GSE232218)
To execute the complete end-to-end differential expression and pathway analysis for Valproic Acid 200 µM:
```bash
py -3.12 scripts/run_vpa_pilot_stage8c.py
```
See the [Reproduction Runbook](REPRODUCTION_RUNBOOK.md) for full instructions and expected statistical outputs.

---

## 4. Documentation Road Map

To explore the project further, consult the following dedicated guides:

| Document | Purpose |
| :--- | :--- |
| **[Project Overview](PROJECT_OVERVIEW.md)** | Scientific foundations, organoid context, and DNT modeling scope |
| **[Current Status](CURRENT_STATUS.md)** | Implemented vs. planned modules, test counts, and release milestones |
| **[Pipeline Map](PIPELINE_MAP.md)** | Visual flowchart of data intake, response contrast, features, and ML |
| **[Data Lineage](DATA_LINEAGE.md)** | Detailed data flow from raw count matrices to model-ready features |
| **[Branch & Release Map](BRANCH_AND_RELEASE_MAP.md)** | Git branch strategy, frozen tags, and release policies |
| **[Reproduction Runbook](REPRODUCTION_RUNBOOK.md)** | Step-by-step instructions to reproduce frozen scientific benchmarks |
| **[Troubleshooting Guide](TROUBLESHOOTING.md)** | Diagnostic matrix covering the 5 data intake failure domains and DE errors |
| **[Release History](RELEASE_HISTORY.md)** | Historical stage-to-scientific-name crosswalk and changelog |
| **[Pipeline Modules](pipeline/data_intake_qc_and_standardization/README.md)** | In-depth module documentation for each pipeline component |
