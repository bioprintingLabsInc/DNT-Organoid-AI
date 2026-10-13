# Cross-Study Gene Response Features

**Historical Identifier**: Stage 8A  
**Scientific Name**: Cross-Study Gene Response Features
**Status**: **STATUS: PLANNED — NO IMPLEMENTATION YET**

---

## 1. Scientific Objective & Governing Design

The purpose of Cross-Study Gene Response Features is to synthesize gene-level molecular perturbation estimates across independent brain-organoid studies without introducing bias or leaking DNT outcome labels.

### Core Methodological Rules:
1. **Within-Study QC Qualification**: Uses strictly QC-qualified within-study treatment-vs-control effects from Differential Molecular Response (`gene_results`).
2. **Preservation of Statistical Metrics**: Preserves effect size ($\log_2\text{FC}$) and its standard error ($\text{SE}$) alongside full experimental context: chemical compound identity, dose/concentration, exposure duration, developmental stage, organoid system, and contrast identity.
3. **No Unweighted Averaging**: Simple unweighted averaging of $\log_2\text{FC}$ values across datasets is prohibited.
4. **No Compound Pooling**: Different chemical compounds are never pooled together.
5. **Experiment Independence**: Multiple doses or exposure timepoints from the same experiment must not be treated as independent studies.
6. **Biologically Comparable Evidence**: Only biologically comparable evidence (comparable developmental stage, organoid context, and exposure window) may be synthesized.
7. **Random-Effects Meta-Analysis**:
   * When multiple independent, biologically comparable studies exist for the same compound regimen, estimates are pooled using **random-effects meta-analysis**.
   * Planned meta-analytic metrics: pooled effect size, pooled standard error/uncertainty, number of supporting independent studies, direction consistency across studies, between-study variance ($\tau^2$), heterogeneity ($I^2$), and leave-one-study-out stability.
8. **Single-Study Preservation**: When only one independent study exists for a compound/regimen (such as the single GSE232218 pilot study), the module retains the study-specific effect size and uncertainty directly without running meta-analysis.

---

## 2. Intended Inputs & Outputs

* **Intended Upstream Input**: Validated `DifferentialExpressionContrastResult` objects from `pipeline/differential_molecular_response/` across one or more independent studies.
* **Intended Output**: Gene-level cross-study feature records preserving:
  * Canonical gene identifier (`canonical_gene_id`, `approved_symbol`);
  * Compound identity and regimen (`agent_identifier`, `concentration`, `exposure_duration`, `developmental_stage`);
  * Single-study estimate (when $N_{\text{studies}} = 1$);
  * Meta-analytic pooled effect, pooled uncertainty, supporting study count, direction consistency, $\tau^2$, and $I^2$ (when $N_{\text{studies}} \ge 2$);
  * Complete source contrast and study provenance IDs.

---

## 3. Implementation Plan

* Implementation will reside in `src/features/cross_study/` and `tests/features/cross_study/`.
* No code paths or algorithms are currently implemented.
