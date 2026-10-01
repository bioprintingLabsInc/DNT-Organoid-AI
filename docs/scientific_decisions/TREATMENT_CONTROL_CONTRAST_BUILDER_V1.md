# Step 7A: Treatment-versus-Matched-Control Contrast Builder v1
## Scientific and Architectural Design Decision

**Status:** APPROVED FOR IMPLEMENTATION
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
   │   ├── Developmental collection age exact matching & exposure separation
   │   ├── Biological replication verification (deflating technical replicates)
   │   ├── Biological & organoid context compatibility verification
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

### 3.1 Explicit Metadata Relationships Only
Under no circumstances does Step 7A infer or guess treatment-control pairings using heuristics (such as substring matching on sample names or assuming unperturbed samples in the same batch are controls).
- Every contrast must originate from an explicitly declared entry in the `treatment_control_relationships` list of the experiment metadata contract.
- Each relationship entry specifies:
  - `treatment_condition_ids` (or `treatment_condition_id`): identifiers for the perturbation conditions.
  - `control_condition_id`: identifier for the baseline/reference control condition.
  - `relationship_type`: semantic role (e.g., `vehicle_control`, `untreated_control`, `solvent_control`).
  - `relationship_id`: unique identifier for the linkage.

### 3.2 One Contrast per Distinct Treatment Condition
When an experimental relationship specifies multiple treatment arms against a shared control (e.g., a concentration series where `treatment_condition_ids = ["VPA_100uM", "VPA_300uM", "VPA_1000uM"]` are matched against `control_condition_id = "Vehicle_DMSO"`):
- The builder constructs **separate, distinct molecular response contrasts**:
  - `VPA_100uM vs. Vehicle_DMSO`
  - `VPA_300uM vs. Vehicle_DMSO`
  - `VPA_1000uM vs. Vehicle_DMSO`
- Under no circumstances are separate treatment conditions pooled into an artificial multi-dose group.
- The shared control samples are correctly associated with each respective contrast, preserving statistical independence across distinct dosing or time levels.

### 3.3 Strict Study and Experiment Boundary Isolation
Cross-study and cross-experiment comparisons are biologically invalid due to un-modeled technical batch variations, differing laboratory protocols, and non-overlapping background noise.
- Every sample in a contrast must belong to the exact same `study_id` and `experiment_id`.
- If a declared relationship spans across distinct studies or experiments, the contrast is immediately marked `ContrastStatus.BLOCKED` with `Severity.ERROR` (`cross_study_comparison_blocked` or `cross_experiment_comparison_blocked`).

### 3.4 Developmental Age and Stage Invariants

Neural organoid development is a highly dynamic temporal process where neurogenesis, gliogenesis, and synaptogenesis progress over days to weeks. Comparing a perturbation at Day 30 against a control at Day 60 would confound chemical perturbation responses with massive developmental maturation gene expression changes.

Step 7A enforces strict temporal invariants:
1. **Exact Matching of Collection Age:**
   - The primary developmental time point for biological comparison is `developmental_age_or_stage_at_collection` (or `developmental_stage_at_collection`).
   - Treatment and control conditions must share identical collection stages (e.g., Day 35 organoids compared against Day 35 controls).
   - Clear mismatches in collection age (e.g., Day 21 vs. Day 42) immediately **BLOCK** the contrast (`developmental_age_mismatch`, `Severity.ERROR`).
2. **Ambiguity and Inconsistency Handling:**
   - If collection age metadata is missing, ambiguous, or discordant among replicates within an arm, the contrast is placed in `ContrastStatus.NEEDS_REVIEW` (`ambiguous_developmental_age`, `Severity.REVIEW`).
   - The platform never guesses, extrapolates, or silently averages developmental stages.
3. **Separation of Exposure Timeline from Collection Stage:**
   - The builder explicitly isolates `exposure_start_developmental_age`, `exposure_duration`, and `washout_duration` from the collection age.
   - For example, an acute exposure (Day 28 to Day 30) and an extended exposure (Day 14 to Day 30) both harvested at Day 30 share a valid collection-age control (Day 30), but maintain distinct exposure metadata blocks for subsequent temporal/pharmacokinetic modeling.
4. **Distinct Collection Stages Produce Distinct Contrasts:**
   - If a multi-stage time-course experiment includes treatments and controls collected at Day 14, Day 28, and Day 56, Step 7A creates distinct contrasts for each collection age (`Treatment_D14 vs Control_D14`, `Treatment_D28 vs Control_D28`, etc.).
   - This design preserves clean pairwise contrasts while enabling future higher-order developmental-age interaction models ($Y \sim \text{Age} + \text{Treatment} + \text{Age}:\text{Treatment}$).

### 3.5 Biological and Organoid Context Compatibility
Perturbations must be compared against controls grown within an identical biological system. The builder enforces:
- **Organoid Context Compatibility:** Matching `organoid_context_id`, organoid type (e.g., cerebral organoid vs. ventral midbrain organoid), and differentiation protocol. Mismatches block contrast formation (`organoid_context_mismatch`, `Severity.ERROR`).
- **Biological Source Consistency:** Matching donor genetic background, cell line (e.g., hiPSC line ID), species (`Homo sapiens`), and sex.

### 3.6 Biological Replication and Inferential Gating

Accurate statistical inference in transcriptomic differential expression depends strictly on biological variance estimation.
1. **Biological Replicate Count ($N_{\text{bio}}$):**
   - The number of biological replicates is defined as the number of **unique `biological_replicate_id` values** present in the sample group.
   - Downstream inferential differential expression requires a minimum of **$N \ge 2$ biological replicates** in both the treatment arm and the control arm.
2. **Technical Replicate Deflation:**
   - Technical replicates (e.g., multiple sequencing libraries or re-sequencing runs of the exact same organoid RNA extraction) share the same `biological_replicate_id`.
   - Technical replicates **never inflate** the biological replicate count ($N_{\text{bio}}$). For instance, four sequencing runs from one organoid constitute $N_{\text{bio}} = 1$.
3. **Inferential Eligibility (`is_inferentially_eligible`):**
   - A contrast is marked `is_inferentially_eligible = True` if and only if:
     $$\text{status} == \text{ContrastStatus.ELIGIBLE} \quad\land\quad N_{\text{bio, treatment}} \ge 2 \quad\land\quad N_{\text{bio, control}} \ge 2$$
   - Contrasts with $N_{\text{bio}} = 1$ in either arm are flagged `under_replicated_condition` (`Severity.REVIEW`), assigned status `NEEDS_REVIEW`, and have `is_inferentially_eligible = False`.
   - Such contrasts are retained for exploratory or descriptive reporting but are strictly barred from un-flagged inferential hypothesis testing.

### 3.7 Explicit Sample Comparison Exception Handling

In complex organoid experiments, specific samples may be flagged by wet-lab technicians or automated QC as compromised (e.g., batch-specific pipetting error, culture vessel contamination, poor viability). The metadata contract allows defining `sample_comparison_exceptions`.

Step 7A processes comparison exceptions with complete auditability:
1. **Disqualification with Valid Replacement:**
   - If an exception flags a control sample as invalid and designates a valid replacement sample from the same cohort and condition, the builder performs the substitution and records an informational finding (`sample_comparison_exception_applied`, `Severity.INFO`).
2. **Disqualification without Replacement:**
   - If an exception disqualifies a sample and no replacement is provided, the builder excludes the sample and flags the contrast as `unresolved_sample_comparison_exception` (`Severity.REVIEW`), placing it into `ContrastStatus.NEEDS_REVIEW`.
   - The builder **never silently drops samples** or silently alters group compositions without a recorded finding.

### 3.8 Input Ordering Invariance and Determinism
Contrast building is fully deterministic.
- The order of sample IDs within contrast groups is canonically sorted.
- The order of contrasts produced in `ResponseContrastDataset` is canonically sorted by `contrast_id`.
- Permuting input sample orders, metadata row orders, or relationship declaration orders produces bit-for-bit identical `ResponseContrastDataset` instances.

---

## 4. Contrast Status and Finding Taxonomy

Each contrast is evaluated against all validation invariants. Findings are recorded as structured, immutable records:
- `name`: Short machine-readable token (e.g., `developmental_age_mismatch`).
- `severity`: One of `Severity.INFO`, `Severity.REVIEW`, or `Severity.ERROR`.
- `message`: Clear, human-readable scientific rationale.
- `metadata`: Key-value context capturing affected samples, IDs, and values.

The final status of the contrast is determined deterministically:
$$\text{Status} = \begin{cases}
\text{BLOCKED}, & \text{if any finding has } \text{severity} == \text{Severity.ERROR} \\
\text{NEEDS\_REVIEW}, & \text{if any finding has } \text{severity} == \text{Severity.REVIEW} \\
\text{ELIGIBLE}, & \text{otherwise}
\end{cases}$$

---

## 5. Architectural Data Models

The implementation in `src/response_builder/models.py` provides the following immutable schemas:

```python
@dataclass(frozen=True)
class Finding:
    name: str
    severity: Severity  # INFO, REVIEW, ERROR
    message: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

@dataclass(frozen=True)
class MolecularResponseContrast:
    contrast_id: str
    study_id: str
    experiment_id: str
    relationship_id: str
    relationship_type: str
    treatment_condition_id: str
    treatment_sample_ids: tuple[str, ...]
    control_condition_id: str
    control_sample_ids: tuple[str, ...]
    developmental_collection_age: str
    organoid_context_id: str
    exposure_metadata: Mapping[str, Any]
    treatment_bio_reps: tuple[str, ...]
    control_bio_reps: tuple[str, ...]
    status: ContrastStatus  # ELIGIBLE, NEEDS_REVIEW, BLOCKED
    findings: tuple[Finding, ...]
    raw_count_matrix_ref: str
    diagnostic_normalized_count_matrix_ref: str

    @property
    def is_inferentially_eligible(self) -> bool:
        return (
            self.status == ContrastStatus.ELIGIBLE
            and self.treatment_replicate_count >= 2
            and self.control_replicate_count >= 2
        )
```

---

## 6. Future Compatibility: Developmental-Age Interaction Models

While Step 7A focuses on pairwise treatment-vs-matched-control contrasts at discrete collection stages, its explicit preservation of:
1. `developmental_age_or_stage_at_collection`,
2. `exposure_start_developmental_age`,
3. `exposure_duration`, and
4. `washout_duration`

guarantees seamless integration into future multi-factor generalized linear models:
$$\log(\mu_{ijk}) = \beta_0 + \beta_{\text{age}} \cdot \text{Age}_j + \beta_{\text{treatment}} \cdot \text{Treatment}_k + \beta_{\text{interaction}} \cdot (\text{Age}_j \times \text{Treatment}_k) + \log(s_j)$$

This ensures the platform can evaluate whether chemical-induced neurodevelopmental toxicity alters the normal developmental trajectory (heterochrony, accelerated senescence, or developmental arrest) rather than merely acting as an acute cytotoxic insult.
