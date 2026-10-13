# Differential Molecular Response

**Historical Identifier**: Step 7B  
**Validation Status**: **VALIDATED** (Commit: `f1a7ad8`)

---

## 1. Scientific & Computational Purpose
Executes negative binomial generalized linear model fitting (Bioconductor DESeq2 Wald test) for an eligible treatment-versus-matched-control contrast. Computes log2 fold changes, Wald statistics, and Benjamini-Hochberg false discovery rates (FDR) while strictly auditing that the entire canonical gene universe is accounted for without silent gene dropping.

### Core Scientific Invariants:
1. **Gene-Universe Conservation**: Every input gene row must be deterministically classified as evaluated canonical or explicitly audited excluded row.
2. **Authoritative DESeq2 Execution**: Fits official DESeq2 dispersion curves and Wald test statistics via a headless R bridge (`DESeq2 1.46.0`).
3. **No Filtering of Non-Significant Genes**: The resulting `gene_results` preserves the entire genome-wide molecular perturbation response, providing unbiased inputs for downstream pathway ranking.

---

## 2. Inputs & Outputs
* **Inputs**: Eligible `MolecularResponseContrast` + count matrix + sample colData.
* **Outputs**: `DifferentialExpressionContrastResult` containing `gene_results: tuple[GeneDifferentialExpressionResult, ...]`, `excluded_genes: tuple[ExcludedGeneAudit, ...]`, and `is_gene_universe_conserved: bool`.

---

## 3. Code, Script & Test Map
* **Primary Code**:
  * [`src/differential_expression/design.py`](../../../../src/differential_expression/design.py): Design formula validation.
  * [`src/differential_expression/models.py`](../../../../src/differential_expression/models.py): Contrast result and gene result dataclasses.
  * [`src/differential_expression/r_bridge.py`](../../../../src/differential_expression/r_bridge.py): Subprocess execution bridge to headless R.
  * [`src/differential_expression/pipeline.py`](../../../../src/differential_expression/pipeline.py): Orchestration engine (`run_contrast_differential_expression`).
* **Headless R Runner**:
  * [`scripts/r/run_deseq2_contrast.R`](../../../../scripts/r/run_deseq2_contrast.R): Official R DESeq2 execution script.
* **Test Suite**:
  * [`tests/differential_expression/test_differential_expression.py`](../../../../tests/differential_expression/test_differential_expression.py)
  * [`tests/differential_expression/test_design.py`](../../../../tests/differential_expression/test_design.py)

---

## 4. Failure Modes & Inspection Protocol
* **Symptom: `is_gene_universe_conserved` is `False`**:
  * *First file to inspect*: [`src/differential_expression/pipeline.py`](../../../../src/differential_expression/pipeline.py) (inspect `run_contrast_differential_expression()`).
  * *Cause*: A gene was dropped during model fitting or count matrix conversion.
* **Symptom: DESeq2 fit error in R bridge**:
  * *First file to inspect*: [`src/differential_expression/r_bridge.py`](../../../../src/differential_expression/r_bridge.py).
  * *Second file to inspect*: [`scripts/r/run_deseq2_contrast.R`](../../../../scripts/r/run_deseq2_contrast.R).
  * *Common cause*: Collinear design matrix or insufficient counts across conditions.
