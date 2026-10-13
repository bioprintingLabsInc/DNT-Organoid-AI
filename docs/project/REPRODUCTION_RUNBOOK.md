# Scientific Reproduction Runbook: GSE232218 VPA Pilot

This runbook provides step-by-step instructions to reproduce the validated Differential Molecular Response and frozen Pathway and Signaling Feature results on public dataset **GSE232218**.

---

## 1. Study and Experimental Context

* **GEO Accession**: `GSE232218`
* **Publication**: Human cortical organoids exposed to chemical compounds during neurodevelopment.
* **Organoid Model**: Cerebral cortex organoids (SOSR-CO), biological source WA09 (H9 hESC).
* **Developmental Age**: Day 6
* **Test Compound**: Valproic Acid (VPA, CID:3121)
* **Concentration**: 200 µM
* **Exposure Duration**: 48 h
* **Vehicle / Control**: Water / vehicle
* **Sample Allocation**: 4 treated biological replicates vs. 4 matched control biological replicates.
  * **Treatment (VPA 200 µM)**: `Sample_3901-AT-3`, `Sample_3901-AT-7`, `Sample_3901-AT-12`, `Sample_3901-AT-17`
  * **Control (Vehicle)**: `Sample_3901-AT-1`, `Sample_3901-AT-5`, `Sample_3901-AT-10`, `Sample_3901-AT-15`

---

## 2. Headless Pilot Execution

Execute the standalone verification runner:
```bash
py -3.12 scripts/run_vpa_pilot_stage8c.py
```

> [!WARNING]
> **Environment & Portability Limitation**:
> The headless pilot script [`scripts/run_vpa_pilot_stage8c.py`](../../scripts/run_vpa_pilot_stage8c.py) currently contains a local-machine path dependency targeting `C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt.gz`.
> * The script is **fully validated and reproducible** in the primary development environment where the downloaded count matrix resides at that location.
> * It is **not yet fully portable to a fresh clone or third-party machine** without placing `GSE232218_gene_expected_count.txt.gz` at that expected path, pending CLI argument parameterization in a future code update.
> * One-command cross-machine reproduction is therefore not claimed until that path parameterization is implemented.

### Expected Output & Acceptance Gates:

1. **Count Matrix & QC**:
   * Total input rows: `60,628`
   * Canonical Ensembl evaluated genes: `60,207`
   * Excluded non-canonical annotation rows: `421`
   * Gene-universe conservation: `True` (evaluated 60,207 + excluded 421 = 60,628 input rows)

2. **Differential Molecular Response (Step 7B)**:
   * Significant upregulated genes ($p_{\text{adj}} < 0.05, \log_2\text{FC} > 0$): **111**
   * Significant downregulated genes ($p_{\text{adj}} < 0.05, \log_2\text{FC} < 0$): **16**

3. **Pathway & Signaling Features (Stage 8C)**:
   * Zero-count genes: `27,980`
   * Genes with valid Wald statistics: `32,227`
   * Cook's distance outliers: `87`
   * Statistically eligible molecular universe: **32,140**
   * Retained non-Entrez valid genes: **7,122** (`future_feature_eligibility=True`, not QC failures)
   * Entrez-mapped eligible genes: `25,018`
   * Duplicate Entrez entries resolved: **2** (Targets: `125775236` and `100874204`)
   * Final ranked genes entering GSEA: **25,016**
   * GO BP pathways evaluated ($15 \le \text{size} \le 500$): `4,893` (1,151 significant at $p_{\text{adj}} < 0.05$)
   * KEGG pathways evaluated ($15 \le \text{size} \le 500$): `348` (129 significant at $p_{\text{adj}} < 0.05$)

4. **Generated Output Artifacts**:
   The script exports the following artifacts to `data/processed/`:
   * `GSE232218_VPA_200uM_statistically_eligible_universe.csv`
   * `GSE232218_VPA_200uM_retained_non_entrez_genes.csv`
   * `GSE232218_VPA_200uM_GO_BP_results.csv`
   * `GSE232218_VPA_200uM_KEGG_results.csv`
   * `GSE232218_VPA_200uM_qc.json`
   * `GSE232218_VPA_200uM_provenance.json`

---

## 3. Interactive Web Interface Reproduction

To reproduce through the browser:
1. Run `py -3.12 -m streamlit run app.py`
2. Open `http://localhost:8501` in your browser.
3. In the sidebar, enter GEO Accession: `GSE232218`.
4. Upload `GSE232218_gene_expected_count.txt.gz`.
5. Verify the 4 treatment samples and 4 control samples are assigned as listed above.
6. Click **Run Differential Expression Analysis**.
7. Confirm the volcano plot renders with 111 up and 16 down DEGs.
