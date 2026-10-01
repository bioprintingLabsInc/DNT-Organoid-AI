# Step 7A: Treatment-versus-Matched-Control Contrast Builder v1
## Scientific and Architectural Design Decision

**Status:** APPROVED FOR IMPLEMENTATION (HARDENED)
**Stage:** Step 7A (Pre-Differential Expression Contrast Formulation)
**Upstream Dependencies:**
- Step 5: Ingestion & Metadata Validation (Metadata Contract)
- Step 6A: Bulk Raw-Count QC (`BulkRawCountQC`)
- Step 6B: Gene Harmonization v1 (Ensembl Release 112 Canonical Human Universe)
- Step 6C: Bulk RNA-seq Normalization v1 (`NormalizedDataset`, `NormalizedCohortResult`, Size Factors $s_j$)

---

## 1. Executive Summary & Scientific Purpose

Step 7A establishes the **Treatment-versus-Matched-Control Contrast Builder v1** for human neural organoid and cell-based transcriptomics within the `DNT-Organoid-AI` platform.

The primary objective of Step 7A is to take the upstream validated datasets—specifically the metadata contract, raw integer count matrix, and normalization cohort outputs from Step 6C—and construct mathematically and biologically sound pairwise comparison contrasts:
$$\text{Contrast}_k = \left(\text{Treatment Condition}_k \text{ vs. Matched Control Condition}_k\right)$$

This stage acts as an unyielding biological and technical gatekeeper. It ensures that every contrast formulated for subsequent differential expression (Step 7B) represents a valid, controlled, un-confounded biological comparison.

### Strict Scope Boundaries

To preserve strict separation of concerns, Step 7A enforces strict scientific and architectural boundaries:
- **No Differential Expression:** Does **not** fit negative binomial generalized linear models (GLMs), estimate dispersions, compute Wald or LRT statistics, or execute Bioconductor `DESeq2`.
- **No Effect-Size Estimation:** Does **not** compute raw or shrunk log2 fold changes ($\log_2\text{FC}$, apeglm, ashr).
- **No Significance Filtering:** Does **not** compute p-values, adjust false discovery rates (FDR / Benjamini-Hochberg), or apply arbitrary significance thresholds.
- **No Gene Filtering:** Does **not** filter or subset genes based on expression levels, variability, fold change, or DNT relevance. The full canonical human gene universe (63,140 Ensembl 112 features) is preserved intact.
- **No DNT Phenotypic or Label Knowledge:** Step 7A is entirely **blind** to downstream DNT classification labels, toxicant categories, reference compounds, known neurodevelopmental genes, or phenotype calls. DNT labels must never influence contrast construction.
- **No Machine Learning or Signature Discovery:** Does **not** perform feature selection, PCA, clustering, pathway enrichment, classifier training, or probability calibration.

---

## 2. Pipeline Positioning and Substrate Immutability

```text
Step 5: Metadata Contract & Schema Validation
   │
Step 6A: Bulk Raw-Count Quality Control (Integer counts, zero-count filtering, sample-level metrics)
   │
Step 6B: Gene Harmonization v1 (Ensembl Release 112 GRCh38.p14 canonical gene universe)
   │
Step 6C: Bulk RNA-seq Normalization v1 (DESeq2 median-of-ratios size factors, cohort verification)
   │
Step 7A: Treatment-versus-Matched-Control Contrast Builder v1 (THIS STAGE)
   │   ├── Condition-level contrast formulation from explicit relationships
   │   ├── Conservation invariant: zero silently dropped treatment conditions
   │   ├── Multi-exposure preservation (ContrastExposure records, deviation isolation)
   │   ├── Developmental age strict matching (prefer governed normalized values)
   │   ├── Organoid context compatibility on scientific attributes (not ID-only)
   │   ├── Biological replication verification (deflating technical replicates)
   │   ├── Explicit sample comparison exception handling
   │   └── Contrast status gating (ELIGIBLE, NEEDS_REVIEW, BLOCKED)
   ▼
Step 7B: Negative Binomial Differential Expression (DESeq2 GLM: raw counts ~ size factors + design)
   │
Step 8: Standardized Molecular Response Feature Engineering
   │
Step 9: Model Training, Conformal Calibration, & DNT Classification
```

### The Immutable Substrate Principle
A foundational tenet of the `DNT-Organoid-AI` architecture is that **raw integer count matrices** generated at Step 6A and harmonized at Step 6B are the **sole primary substrate** for all downstream statistical modeling.

1. **Raw Counts for Differential Expression:** Downstream negative-binomial modeling (DESeq2) expects integer counts as inputs and accounts for sequencing depth differences via the size factors ($s_j$) estimated in Step 6C.
2. **Diagnostic Normalized Counts for Technical Inspection Only:** The normalized count matrix calculated in Step 6C ($K_{ij} / s_j$) is strictly diagnostic. It is utilized exclusively for sample visualization, outlier auditing, and technical QC checks. It is **never** supplied as input counts to downstream differential expression models.

---

## 3. Core Architectural and Scientific Invariants

### 3.1 No Silently Dropped Treatment Conditions (Conservation Invariant)
Every treatment condition represented in an eligible normalized cohort must have an explicit, auditable disposition:
- **`ELIGIBLE`** contrast
- **`NEEDS_REVIEW`** contrast/disposition
- **`BLOCKED`** contrast/disposition

A treatment condition must **never disappear** merely because:
- `treatment_control_relationship` is missing;
- `matched_control_condition_ids` is empty;
- linked control samples are missing;
- the relationship is malformed; or
- a comparison exception prevents normal matching.

If a treatment condition exists in a normalized cohort but lacks a valid declared relationship:
- It produces an explicit contrast with `status = ContrastStatus.BLOCKED`.
- It records `rule_id = "missing_explicit_treatment_control_relationship"`.
- It preserves `study_id`, `experiment_id`, `cohort_id`, `treatment_condition_id`, and `treatment_sample_ids`.

**Conservation Invariant:**
$$\text{Total Treatment Conditions Encountered} = \text{Eligible} + \text{Needs Review} + \text{Blocked Dispositions}$$
No treatment condition may be lost from the audit trail.

### 3.2 One Contrast per Distinct Treatment Condition
When an experimental relationship specifies multiple treatment arms against a shared control (e.g., `treatment_condition_ids = ["VPA_100uM", "VPA_300uM", "VPA_1000uM"]` matched against `control_condition_id = "Vehicle_DMSO"`):
- The builder constructs **separate, distinct molecular response contrasts**:
  - `VPA_100uM vs. Vehicle_DMSO`
  - `VPA_300uM vs. Vehicle_DMSO`
  - `VPA_1000uM vs. Vehicle_DMSO`
- Separate treatment conditions are never pooled into an artificial multi-dose arm.
- Shared control samples are associated with each respective contrast, preserving statistical independence.

### 3.3 Preservation of Multiple Exposure Records (`ContrastExposure`)
Complex perturbations frequently involve multi-agent or factorial regimens (e.g., drug + cytokine, co-exposure, factorial challenge).
- The contrast preserves **all condition-planned exposure records** in `treatment_exposures: tuple[ContrastExposure, ...]`.
- Exposures are sorted deterministically by `(exposure_id, agent_name, agent_identifier)`.
- Backward-compatible singular fields (`agent_name`, `agent_identifier`, `vehicle`, etc.) are populated **only** when there is exactly one applicable exposure (`len(treatment_exposures) == 1`). If multiple exposures exist, singular fields are set to `None` so they **never hide a multi-agent condition**.
- Sample-level exposure deviations (records with `sample_id`) are preserved separately in `sample_level_exposure_deviations` with an informational finding and **never overwrite** condition-planned exposures.

### 3.4 Strict Study and Experiment Boundary Isolation
Cross-study and cross-experiment comparisons are biologically invalid due to un-modeled technical batch variations, differing laboratory protocols, and non-overlapping background noise.
- Every sample in a contrast must belong to the exact same `study_id` and `experiment_id`.
- Conflicting study or experiment IDs immediately trigger `ContrastStatus.BLOCKED` (`cross_study_comparison_blocked` or `cross_experiment_comparison_blocked`).

### 3.5 Developmental Age and Stage Invariants
Neural organoid development is a highly dynamic temporal process where neurogenesis, gliogenesis, and synaptogenesis progress over days to weeks. Comparing a perturbation at Day 30 against a control at Day 60 would confound chemical perturbation responses with massive developmental maturation gene expression changes.

Step 7A enforces strict temporal invariants:
1. **Governed Normalized Age Preferred:**
   - When a canonical `value_record` provides a governed `normalized_value` (e.g. `day_30`), that normalized representation is preferred over raw text strings.
   - Scientifically equivalent representations (such as `"Day 30"`, `"30 days"`, `30`) mapped to `day_30` are accepted as matching.
2. **True Known Age Mismatch Blocks:**
   - Clear mismatches in collection age (e.g., Day 30 vs. Day 60) immediately **BLOCK** the contrast (`developmental_age_mismatch`, `Severity.ERROR`).
3. **Unresolved Text Triggers Review (Never Guessed):**
   - If an age string cannot be standardized to a known temporal unit and lacks a governed `normalized_value` (e.g. `"mid-maturation"` vs. `"stage 2"`), the builder **never invents a normalized age**.
   - It records `rule_id = "unresolved_developmental_age_equivalence"` with `Severity.REVIEW`, placing the contrast into `ContrastStatus.NEEDS_REVIEW`.
4. **Separation of Exposure Timeline from Collection Stage:**
   - `developmental_age_or_stage_at_exposure`, `exposure_start_time_or_stage`, `exposure_duration`, and `washout_or_recovery_duration` are tracked separately from `treatment_collection_age_or_stage` and `control_collection_age_or_stage`.

### 3.6 Organoid Context Compatibility (Scientific Attributes over ID Identity)
Perturbations must be compared against controls grown within a biologically compatible organoid system.
- Compatibility is determined by **governed scientific attributes**, not identifier equality alone:
  - `organoid_type` (e.g. cerebral organoid vs hepatic organoid)
  - `brain_region_or_model_identity` (e.g. cerebral cortex vs ventral midbrain)
  - `differentiation_protocol_identifier`
  - `differentiation_protocol_version`
- **Compatible Contexts with Different IDs:** If treatment and control have different `organoid_context_id` values (e.g. `"ctx_treat_1"` and `"ctx_ctrl_1"`), but their governed attributes match, the comparison **remains eligible**.
- **Context Mismatch:** Conflicting `organoid_type` or `brain_region_or_model_identity` immediately **BLOCKS** the contrast (`organoid_context_mismatch`).
- **Unresolved Context:** Missing context records or blank attributes trigger `ContrastStatus.NEEDS_REVIEW` (`unresolved_organoid_context_compatibility`).

### 3.7 Biological Source Handling
- Unpaired multi-donor or multi-cell-line designs within the same study and experiment are biologically valid and remain `ELIGIBLE`.
- All `biological_source_ids` and `biological_replicate_ids` are preserved on the contrast.
- Conflicting species (e.g. human vs mouse) immediately **BLOCKS** the contrast (`cross_species_comparison_blocked`).

### 3.8 Biological Replication and Inferential Gating
Accurate statistical inference in transcriptomic differential expression depends strictly on biological variance estimation.
1. **Biological Replicate Count ($N_{\text{bio}}$):**
   - Defined as the number of **unique `biological_replicate_id` values** present in the sample group.
   - Downstream inferential differential expression requires a minimum of **$N \ge 2$ biological replicates** in both the treatment arm and the control arm.
2. **Technical Replicate Deflation:**
   - Technical replicates sharing the same `biological_replicate_id` **never inflate** $N_{\text{bio}}$.
3. **Inferential Eligibility (`is_inferentially_eligible`):**
   - A contrast is marked `is_inferentially_eligible = True` if and only if:
     $$\text{status} == \text{ContrastStatus.ELIGIBLE} \quad\land\quad N_{\text{bio, treatment}} \ge 2 \quad\land\quad N_{\text{bio, control}} \ge 2$$
   - Contrasts with $N_{\text{bio}} = 1$ in either arm are flagged `insufficient_biological_replication` (`Severity.REVIEW`), assigned status `NEEDS_REVIEW`, and have `is_inferentially_eligible = False`.

### 3.9 Explicit Sample Comparison Exception Handling
- Disqualification with valid replacement: replaces affected control samples, records `sample_comparison_exception_applied` (`Severity.INFO`), remains eligible if other criteria met.
- Disqualification without replacement: flags `unresolved_sample_comparison_exception` (`Severity.REVIEW`), moves to `ContrastStatus.NEEDS_REVIEW`.
- Samples are **never silently dropped**.

### 3.10 Input Ordering Invariance and Determinism
- Contrast sample IDs and biological replicate IDs are canonically sorted.
- Contrast exposures are sorted by `(exposure_id, agent_name, agent_identifier)`.
- Contrasts in `ResponseContrastDataset` are sorted by `(study_id, experiment_id, cohort_id, treatment_condition_id, contrast_id)`.
- Permuting input orders produces identical contrast objects and hashes.

---

## 4. Architectural Data Models

The implementation in `src/response_builder/models.py` provides:

```python
@dataclass(frozen=True, slots=True)
class ContrastExposure:
    exposure_id: str
    agent_name: str | None = None
    agent_identifier: str | None = None
    vehicle: str | None = None
    concentration_or_dose: float | str | None = None
    concentration_or_dose_unit: str | None = None
    exposure_start_time_or_stage: str | None = None
    developmental_age_or_stage_at_exposure: str | None = None
    exposure_duration: float | str | None = None
    exposure_duration_unit: str | None = None
    washout_or_recovery_duration: float | str | None = None
    washout_or_recovery_duration_unit: str | None = None

@dataclass(frozen=True, slots=True)
class MolecularResponseContrast:
    contrast_id: str
    study_id: str
    experiment_id: str
    cohort_id: str
    relationship_id: str
    treatment_condition_id: str
    matched_control_condition_ids: tuple[str, ...]
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]
    treatment_biological_replicate_ids: tuple[str, ...]
    control_biological_replicate_ids: tuple[str, ...]
    treatment_exposures: tuple[ContrastExposure, ...]
    sample_level_exposure_deviations: tuple[dict[str, Any], ...]
    agent_name: str | None
    agent_identifier: str | None
    vehicle: str | None
    concentration_or_dose: float | str | None
    concentration_or_dose_unit: str | None
    exposure_start_time_or_stage: str | None
    developmental_age_or_stage_at_exposure: str | None
    exposure_duration: float | str | None
    exposure_duration_unit: str | None
    washout_or_recovery_duration: float | str | None
    washout_or_recovery_duration_unit: str | None
    treatment_collection_age_or_stage: str | None
    control_collection_age_or_stage: str | None
    treatment_collection_age_normalized: str | None
    control_collection_age_normalized: str | None
    biological_source_ids: tuple[str, ...]
    organoid_context_ids: tuple[str, ...]
    sample_comparison_exception_ids: tuple[str, ...]
    status: ContrastStatus
    findings: tuple[Finding, ...]

    @property
    def is_inferentially_eligible(self) -> bool:
        return (
            self.status == ContrastStatus.ELIGIBLE
            and self.treatment_replicate_count >= 2
            and self.control_replicate_count >= 2
        )
```

---

## 5. Future Compatibility: Developmental-Age Interaction Models

While Step 7A focuses on pairwise treatment-vs-matched-control contrasts at discrete collection stages, its explicit preservation of:
1. `treatment_collection_age_or_stage` and `treatment_collection_age_normalized`,
2. `developmental_age_or_stage_at_exposure`,
3. `exposure_duration`, and
4. `washout_or_recovery_duration`

guarantees seamless integration into future multi-factor generalized linear models:
$$\log(\mu_{ijk}) = \beta_0 + \beta_{\text{age}} \cdot \text{Age}_j + \beta_{\text{treatment}} \cdot \text{Treatment}_k + \beta_{\text{interaction}} \cdot (\text{Age}_j \times \text{Treatment}_k) + \log(s_j)$$

This ensures the platform can evaluate whether chemical-induced neurodevelopmental toxicity alters the normal developmental trajectory (heterochrony, accelerated senescence, or developmental arrest) rather than merely acting as an acute cytotoxic insult.
