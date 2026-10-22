# Cross-Study Gene Response Features

**Historical Identifier**: Stage 8A
**Scientific Name**: Cross-Study Gene Response Features
**Current Milestone**: **Stage 8A.1 — Condition-Level Gene Response Registry (Implemented & Verified on GSE232218 Pilot)**
**Current Status**: **STAGE_8A1_VALIDATED_GSE232218_PILOT**
**Future Milestone**: **Stage 8A.2 — Multi-Study Synthesis & Comparability Evaluation (Planned — Not Implemented)**

---

## 1. Scientific Purpose & Governance Boundary

The purpose of Cross-Study Gene Response Features is to synthesize gene-level molecular perturbation estimates across independent brain-organoid studies without introducing bias or leaking DNT outcome labels.

> [!CRITICAL]
> **Foundational Governance Principle**:
> **Stage 8A.1 creates a verified condition-level gene-response registry on a single-study pilot (GSE232218). It does not by itself establish cross-study reproducibility.**
> Genuine cross-study synthesis is **NOT** yet implemented. Cross-study synthesis (Stage 8A.2) will occur only after a second independent dataset is processed through Stages 1–7 and approved under formal biological and statistical comparability criteria.

### Fundamental Unit of Stage 8A.1
The fundamental scientific unit of Stage 8A.1 is:
$$\mathbf{one\ gene \times one\ treatment\text{-}control\ contrast}$$

Every gene evaluated by Stage 7B in a treatment-versus-matched-control contrast is preserved with its exact differential statistics and complete experimental, chemical, and biological context. Stage 8A.1 strictly preserves condition-level identity rather than collapsing, averaging, or filtering records.

---

## 2. Non-Negotiable Scientific Rules

1. **Exact Numerical & Context Preservation**:
   * Preserves effect size ($\log_2\text{FC}$), standard error ($\text{SE}$), Wald statistic, $p$-value, and FDR/adjusted $p$-value exactly as calculated by Stage 7.
   * Stage 8A never recalculates differential expression statistics.
2. **No Unweighted Averaging**:
   * Simple unweighted averaging of $\log_2\text{FC}$ values across conditions or datasets is strictly prohibited.
3. **No Compound Pooling**:
   * Different chemical compounds are never pooled together.
4. **Study Independence Definition**:
   * Multiple concentrations, doses, timepoints, or contrasts from the same experiment or study accession (e.g. GSE232218 VPA 200 µM and VPA 400 µM) share the same `independent_study_id`.
   * They **must not** be treated as independent studies and do not inflate the independent study count.
5. **Single-Study Status (`SINGLE_STUDY_ONLY`)**:
   * Because the initial foundation uses GSE232218, all gene/compound records are assigned `CrossStudyStatus.SINGLE_STUDY_ONLY`.
   * The pipeline strictly prohibits generating pooled effects, between-study heterogeneity statistics ($\tau^2$, $I^2$), or reproducibility claims from a single study.
6. **No Premature Meta-Analysis Commitment / Method-Neutral Future Contract**:
   * Stage 8A.2 will first determine whether multiple independent studies for the same compound and gene are biologically and statistically comparable. If synthesis is scientifically justified, an appropriate statistical synthesis method will be selected and prespecified at that stage based on the observed study structure, dependence, heterogeneity, number of studies, and experimental context. Stage 8A.1 does not lock fixed-effects, random-effects, DerSimonian-Laird, REML, Hartung-Knapp, or any other estimator.
7. **No Arbitrary Comparability Thresholds**:
   * The pipeline does not invent arbitrary thresholds (e.g. "doses within 2× are comparable" or "stages within X days are comparable").
8. **No Significance or Magnitude Filtering**:
   * Non-significant genes (FDR $> 0.05$) and small-effect genes ($|\log_2\text{FC}| < 0.2$) are **not** discarded. All Stage 7 evaluated genes remain present.
9. **No Entrez ID Requirement**:
   * Canonical Ensembl gene IDs from Ensembl Release 112 remain the primary identifiers. Valid non-Entrez genes (e.g. lncRNAs) are fully preserved.
10. **Zero DNT Label Leakage**:
    * Operates with zero dependency on DNT-positive, DNT-negative, or regulatory toxicology outcome labels.

---

## 3. Data Contracts & Output Schema

### Upstream Input Contract
Consumes validated `DifferentialExpressionContrastResult` (or `DifferentialExpressionDataset`) from `src/differential_expression/`.

### Output Registry Schema (`cross_study_gene_response_long`)
Each row represents one gene evaluated in one treatment-control contrast:

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `study_id` | `str` | Study accession (e.g. `GSE232218`) |
| `study_family_id` | `str` | Study family identifier for multi-experiment series |
| `experiment_id` | `str` | Specific experiment identifier within the study |
| `contrast_id` | `str` | Unique contrast identifier |
| `compound_id` | `str` | Normalized chemical identifier (e.g. `CID:3121`) |
| `compound_name` | `str` | Chemical agent name (e.g. `Valproic Acid`) |
| `dose_value` | `float \| None` | Numeric concentration/dose value (e.g. `200.0`) |
| `dose_unit` | `str \| None` | Concentration unit (e.g. `uM`) |
| `dose_original` | `str \| None` | Verbatim dose string (e.g. `200.0 uM`) |
| `developmental_stage` | `str \| None` | Organoid age/stage (e.g. `Day 6`) |
| `exposure_duration` | `float \| None` | Exposure duration value (e.g. `48.0`) |
| `model_type` | `str \| None` | Organoid model system (e.g. `Cerebral Cortex Organoid`) |
| `canonical_ensembl_gene_id`| `str` | Primary canonical Ensembl ID (e.g. `ENSG00000000003`) |
| `gene_symbol` | `str \| None` | Approved gene symbol (e.g. `TSPAN6`) |
| `base_mean` | `float` | DESeq2 normalized baseMean across contrast samples |
| `log2_fold_change` | `float \| None` | Preserved log2 fold change from Stage 7 |
| `log2_fold_change_se` | `float \| None` | Preserved log2FC standard error from Stage 7 |
| `wald_statistic` | `float \| None` | Preserved Wald statistic from Stage 7 |
| `p_value` | `float \| None` | Preserved two-tailed Wald p-value from Stage 7 |
| `adjusted_p_value` | `float \| None` | Preserved Benjamini-Hochberg FDR from Stage 7 |
| `response_direction` | `str` | Descriptive direction: `UP`, `DOWN`, `NO_DIRECTION`, `UNDETERMINED` |
| `estimate_status` | `str` | Estimation state: `OK`, `ZERO_COUNTS`, `COOKS_OUTLIER`, `FIT_ERROR` |
| `missing_reason` | `str \| None` | Explicit reason if statistical values are missing |
| `independent_study_id` | `str` | Normalized independent study identifier (e.g. `GSE232218`) |
| `cross_study_status` | `str` | Synthesis status: `SINGLE_STUDY_ONLY` |
| `provenance` | `str` | Serialized JSON tracking upstream hashes and versions |

---

## 4. Reconciled Accounting on GSE232218 Pilot

Audited directly on public **GSE232218** (Valproic Acid 200 µM in Day 6 human cortical organoids, 4 treated vs. 4 matched vehicle controls):

```text
Stage 7 Canonical Evaluated Genes                    60,207
└── Stage 8A.1 Condition Records Produced            60,207 (100% 1:1 Preservation)
    |-- Valid Statistically Tested Estimates:        32,140 (EstimateStatus.OK)
    |-- Zero-Count Genes (baseMean == 0):            27,980 (EstimateStatus.ZERO_COUNTS)
    \-- Cook's Distance Outliers:                        87 (EstimateStatus.COOKS_OUTLIER)

Direction Breakdown (Descriptive Only):
    |-- UP (log2FC > 0):                             16,284
    |-- DOWN (log2FC < 0):                           15,943
    |-- NO_DIRECTION (log2FC == 0):                       0
    \-- UNDETERMINED (log2FC is None / NaN):         27,980 (Zero-count genes)

Multi-Study Synthesis Boundary:
    \-- SINGLE_STUDY_ONLY records:                   60,207 (100% of records)
```

### Verification Findings
* **Numerical Preservation**: 0 mismatches across all 60,207 genes for $\log_2\text{FC}$, $\text{SE}$, Wald, $p$, and FDR.
* **Gene Universe Conservation**: 100% of canonical genes from Stage 7 are preserved in the registry.
* **Missing-Data Integrity**: 27,980 zero-count genes retain `None`/`NaN` for $\log_2\text{FC}$ and $\text{SE}$ with explicit `missing_reason` rather than being fabricated as `0.0`.
* **Zero DNT Leakage**: 0 DNT labels or regulatory classifications accessed or present.

---

## 5. Architecture & Implementation Code Paths

* **Public API**: [`src/features/cross_study/__init__.py`](../../../../src/features/cross_study/__init__.py)
* **Data Models & Enums**: [`src/features/cross_study/models.py`](../../../../src/features/cross_study/models.py)
* **Record Builder**: [`src/features/cross_study/builder.py`](../../../../src/features/cross_study/builder.py)
* **Registry Construction**: [`src/features/cross_study/registry.py`](../../../../src/features/cross_study/registry.py)
* **Quality Control Auditor**: [`src/features/cross_study/qc.py`](../../../../src/features/cross_study/qc.py)
* **Serialization Services**: [`src/features/cross_study/serialization.py`](../../../../src/features/cross_study/serialization.py)
* **Automated Test Suite**: [`tests/features/cross_study/`](../../../../tests/features/cross_study/) (25 passing unit and integration tests)
* **Pilot Execution Runner**: [`scripts/run_vpa_pilot_stage8a.py`](../../../../scripts/run_vpa_pilot_stage8a.py)

---

## 6. Future Stage 8A.2 Comparability Boundary (Planned — Not Implemented)

Stage 8A.2 will be designed when a second independent dataset (e.g. an independent organoid study testing Valproic Acid or another compound) completes Stages 1–7.

Stage 8A.2 will first determine whether multiple independent studies for the same compound and gene are biologically and statistically comparable. If synthesis is scientifically justified, an appropriate statistical synthesis method will be selected and prespecified at that stage based on the observed study structure, dependence, heterogeneity, number of studies, and experimental context. Stage 8A.1 does not lock fixed-effects, random-effects, DerSimonian-Laird, REML, Hartung-Knapp, or any other estimator.

### Retained Core Governance Invariants for Stage 8A.2:
* **Separation of Non-Comparable Conditions**: Non-comparable conditions (discrepant developmental age, exposure duration, model system) remain strictly separate and are never forcibly synthesized.
* **Permanent Preservation of Condition Records**: Underlying Stage 8A.1 condition-level records remain permanently preserved and are never collapsed or discarded.
* **No Automatic Synthesis**: No synthesis occurs merely because two studies test the same nominal chemical compound.
* **No Premature Implementation**: Genuine cross-study synthesis is not implemented in Stage 8A.1; GSE232218 is evaluated strictly as a single-study foundation.

Planned multi-study evaluation states:
* `SINGLE_STUDY_ONLY`: Retained as single-study condition record (current verified status).
* `MULTI_STUDY_CONTEXT_UNDER_REVIEW`: Multiple independent studies exist for compound; biological comparability under evaluation.
* `NONCOMPARABLE_CONTEXT`: Discrepant developmental age, exposure window, or organoid platform prevents pooling.
* `COMPARABLE_INDEPENDENT_STUDIES`: Rigorously verified biologically and statistically comparable independent studies eligible for prespecified synthesis evaluation.
* `INSUFFICIENT_METADATA`: Incomplete regimen parameters prevent comparability determination.
