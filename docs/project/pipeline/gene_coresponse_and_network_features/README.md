# Stage 8B: Gene Co-Response and Biological Network Features

**Historical Identifier**: Stage 8B  
**Scientific Name**: Gene Co-Response and Biological Network Features  
**Status**: **IMPLEMENTED & SCIENTIFICALLY GOVERNED**  
**Module Location**: `src/features/coresponse/`  
**Test Suite**: `tests/features/coresponse/`  

---

## 1. Scientific Objective & Scope

Stage 8B organizes, analyzes, and quantifies biological relationships among genes participating in human neurodevelopmental responses to chemical exposure. 

Unlike naive network pipelines that treat all molecular interactions identically, Stage 8B enforces a strict four-layer ontological distinction among evidence types:
1. **Established Biological Reference Networks** (`ESTABLISHED_REFERENCE_INTERACTION`): Pre-existing, external human biological interactomes (e.g. STRING, OmniPath, Reactome PPIs) independent of empirical chemical perturbation cohorts.
2. **Observed Within-Study Gene Co-Expression** (`WITHIN_STUDY_CO_EXPRESSION`): Pairwise continuous expression correlations across biological samples within a single experimental study or cohort.
3. **Cross-Condition Gene Co-Response** (`CROSS_CONDITION_CO_RESPONSE`): Concordant differential expression profiles across multiple independent chemical perturbation contrasts.
4. **Data-Adaptive Prediction Feature Networks** (`DATA_ADAPTIVE_PREDICTION_FEATURE`): Machine-learned or empirical networks explicitly restricted to training partitions to prevent data leakage into AI/ML evaluation.

These four layers are never conflated, combined, or treated as interchangeable evidence.

---

## 2. Statistical Gating & Independence Safeguards

To prevent the generation of spurious correlations or false network hubs from under-powered cohorts, Stage 8B enforces strict sample-size and diversity prerequisites:

### A. Within-Study Co-Expression Prerequisites
- **Independent Biological Samples**: Requires $N_{\text{bio}} \ge 10$ biological replicate samples within the study.
- **Technical Replicate Protection**: Replicates sharing the same biological donor or culture preparation cannot artificially inflate sample size; technical replicates are strictly aggregated or collapsed prior to covariance estimation.
- **Under-Powered Datasets**: Studies with $N_{\text{bio}} < 10$ (such as GSE232218 with $N=4$ per group, or GSE187006 with $N=3$ per group) honestly yield `INSUFFICIENT_SAMPLES` status with zero manufactured edges.

### B. Cross-Condition Co-Response Prerequisites
- **Multi-Contrast & Multi-Chemical Diversity**: Requires $\ge 8$ perturbation contrasts across $\ge 4$ distinct chemical compounds.
- **Independence Enforcement**: Multiple concentrations (e.g. dose-response curves) or timepoints of the same chemical compound do **not** constitute independent chemical evidence.
- **Single-Chemical Cohorts**: Registries with $< 4$ distinct chemical agents (such as GSE232218 and GSE187006 which only test Valproic Acid) honestly emit `INSUFFICIENT_CHEMICALS` status with zero manufactured edges.

---

## 3. AI/ML Partitioning & Data Leakage Prevention

Empirical networks intended as prediction features must never be influenced by held-out test or validation data:
- **Training-Only Partitioning**: `filter_registry_by_training_partition` isolates records exclusively to approved training chemicals and study families.
- **Audit Verification**: `verify_partition_isolation` systematically inspects biological networks to guarantee zero presence of held-out compounds.
- **Supervisory Outcome Quarantine**: Terms indicative of downstream labels (e.g., `dnt_positive`, `dnt_negative`, `gold_standard_positive`) are strictly banned from feature networks.

---

## 4. Network Topology & Scientific Interpretation Guardrails

Topological metrics are extracted using NetworkX:
- **Metrics Computed**: Degree, clustering coefficient, neighborhood size, isolation indicator, and empirical degree percentile rank.
- **Literature Hub Degree Bias**: Well-studied genes (e.g. *TP53*, *EGFR*, *AKT1*) have high degrees in reference databases due to publication volume rather than intrinsic neurodevelopmental specificity. Every topology record attaches the mandatory disclaimer:
  > *"Gene degree reflects reference literature investigation density as well as biological connectivity. High degree does NOT establish causal criticality, direct chemical binding, or validated DNT biomarker status."*
- **Causality & Binding Disclaimers**: Network edges and subnetwork extractions explicitly disclaim causal regulatory control, direct compound target binding, or biomarker validity in the absence of targeted biophysical validation.

---

## 5. Canonical Governed Schemas

### A. 26-Column Canonical Network Edge Schema (`CANONICAL_NETWORK_EDGE_COLUMNS`)
| Column | Description |
|---|---|
| `edge_id` | Deterministic SHA-256 hash of sorted gene pair, relationship type, and evidence source. |
| `gene_a_canonical_id` | Ensembl gene ID (lexicographically first). |
| `gene_b_canonical_id` | Ensembl gene ID (lexicographically second). |
| `gene_a_symbol` | Human HGNC gene symbol for Gene A. |
| `gene_b_symbol` | Human HGNC gene symbol for Gene B. |
| `relationship_type` | One of the 4 non-interchangeable evidence layers. |
| `interaction_nature` | Nature of interaction (e.g. `PHYSICAL_INTERACTION`, `STATISTICAL_CO_RESPONSE`). |
| `evidence_source` | Authority or pipeline module (e.g. `STRING`, `Stage8B_Coresponse`). |
| `evidence_version` | Upstream reference database or algorithm version. |
| `evidence_category` | Derivation category (e.g. `DATABASE_CURATED`, `EMPIRICAL_CORESPONSE`). |
| `confidence_score` | Normalized confidence score in $[0.0, 1.0]$. |
| `statistical_metric_name` | Name of statistic (e.g. `spearman_rho`, `pearson_r`). |
| `statistical_estimate` | Value of correlation or association estimate. |
| `p_value` | Parametric two-tailed $p$-value from Student's $t$-distribution. |
| `adjusted_p_value` | Benjamini-Hochberg False Discovery Rate (FDR). |
| `sample_size` | Number of independent samples used for within-study co-expression. |
| `condition_count` | Number of distinct perturbation contrasts used for co-response. |
| `study_context` | Tuple/list of contributing study identifiers. |
| `compound_context` | Tuple/list of contributing compound identifiers. |
| `partition_assignment` | Safety partition (e.g. `UNRESTRICTED_REFERENCE`, `TRAINING_ONLY`). |
| `qc_status` | Statistical and data-integrity quality control status. |
| `is_direct_binding` | Boolean flag (True only for verified physical protein binding). |
| `is_chemical_specific` | Boolean flag (True for empirical perturbation-derived networks). |
| `causality_claim` | Strict non-causal characterization. |
| `scientific_limitations` | Methodological and statistical boundaries. |
| `provenance_notes` | Complete execution metadata and input hashes. |

### B. 9-Column Canonical Topology Metrics Schema (`CANONICAL_TOPOLOGY_METRICS_COLUMNS`)
1. `canonical_ensembl_gene_id`
2. `gene_symbol`
3. `network_id`
4. `degree`
5. `clustering_coefficient`
6. `neighborhood_size`
7. `is_isolated`
8. `degree_percentile`
9. `reference_degree_bias_note`

---

## 6. Implementation Architecture

```
src/features/coresponse/
├── __init__.py                # Governed public exports
├── models.py                  # Dataclasses and strictly typed Enums
├── reference_network.py       # External reference network ingestion & multi-version support
├── within_study.py            # Within-study co-expression & sample-size gating (N >= 10)
├── cross_condition.py         # Cross-condition co-response & chemical diversity gating
├── topology.py                # NetworkX topology metrics & degree bias governance
├── partitioning.py            # AI/ML partition safety & leakage auditing
├── qc.py                      # Comprehensive network quality control auditor
├── serialization.py           # Governed DataFrame & CSV export functions
└── runner.py                  # High-level pipeline execution orchestrator
```
