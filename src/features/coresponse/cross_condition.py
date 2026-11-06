"""Stage 8B: Cross-Condition Gene Co-Response Analysis and Independence Gating.

Evaluates coordinated gene transcriptional shifts across qualified Stage 8A perturbation
contrasts while strictly enforcing:
1. Diversity & Independence Gating: Requires at least `min_conditions` (default 8) across
   `min_chemicals` (default 4) independent chemical treatments.
   Under-conditioned cohorts (such as single-chemical registries with few contrasts)
   honestly emit INSUFFICIENT_EVIDENCE / INSUFFICIENT_CONDITIONS rather than manufacturing
   artificial correlations or network edges.
2. Metadata Preservation: Retains chemical, dose, duration, model, and study provenance.
3. Separation of Descriptive Concordance from Statistically Supported Co-Response.
4. Non-Interactivity Safeguard: Genes sharing response direction are NOT claimed to physically interact.
"""

from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
import math
from typing import Any
import numpy as np
import pandas as pd
from scipy import stats

from src.features.coresponse.models import (
    BiologicalNetwork,
    CoresponseEvaluationSummary,
    CoresponseQCStatus,
    EvidenceCategory,
    GeneNetworkEdge,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
    derive_edge_id,
)
from src.features.coresponse.within_study import compute_benjamini_hochberg_fdr
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    EstimateStatus,
    ResponseDirection,
)


def validate_cross_condition_eligibility(
    registry: ConditionGeneResponseRegistry,
    min_conditions: int = 8,
    min_chemicals: int = 4,
    training_chemicals: set[str] | None = None,
) -> tuple[bool, int, int, str]:
    """Validate whether an experimental registry contains sufficient independent chemical conditions.

    Args:
        registry: Validated Stage 8A ConditionGeneResponseRegistry.
        min_conditions: Minimum distinct perturbation contrasts required (default 8).
        min_chemicals: Minimum distinct chemical compounds required (default 4).
        training_chemicals: Optional partition filter to restrict eligibility check to training split.

    Returns:
        tuple of (is_eligible: bool, condition_count: int, chemical_count: int, reason: str)
    """
    records = registry.records
    if training_chemicals is not None:
        records = tuple(r for r in records if r.compound_id in training_chemicals or r.compound_name in training_chemicals)

    if not records:
        return False, 0, 0, "Registry contains no valid records matching criteria."

    distinct_contrasts = {r.contrast_id for r in records}
    distinct_chemicals = {r.compound_id for r in records}

    n_cond = len(distinct_contrasts)
    n_chem = len(distinct_chemicals)

    if n_cond < min_conditions:
        return (
            False,
            n_cond,
            n_chem,
            f"Insufficient independent conditions: found {n_cond} contrast(s), but minimum required "
            f"is min_conditions={min_conditions}. Computing cross-condition covariance across few "
            f"perturbations is statistically invalid."
        )

    if n_chem < min_chemicals:
        return (
            False,
            n_cond,
            n_chem,
            f"Insufficient chemical diversity: found {n_chem} chemical(s), but minimum required "
            f"is min_chemicals={min_chemicals}. Multiple conditions of the same compound or "
            f"related cell lines do not constitute independent evidence for general co-response networks."
        )

    return True, n_cond, n_chem, "Statistically eligible for cross-condition co-response inference."


def analyze_cross_condition_coresponse(
    registry: ConditionGeneResponseRegistry,
    network_id: str | None = None,
    network_name: str | None = None,
    min_conditions: int = 8,
    min_chemicals: int = 4,
    method: str = "spearman",
    metric_column: str = "log2_fold_change",  # or "wald_statistic"
    min_abs_correlation: float = 0.7,
    max_fdr: float = 0.05,
    max_evaluated_genes: int | None = 2000,
    training_chemicals: set[str] | None = None,
    partition_assignment: PartitionAssignment = PartitionAssignment.DESCRIPTIVE_EXPLORATION,
) -> tuple[BiologicalNetwork, CoresponseEvaluationSummary]:
    """Execute cross-condition gene co-response analysis with strict independence gating.

    If the registry lacks adequate condition or chemical diversity (e.g. a single contrast, or
    multiple contrasts of only one chemical), halts inference honestly and returns an empty network with
    INSUFFICIENT_EVIDENCE / INSUFFICIENT_CONDITIONS status.

    Args:
        registry: Validated Stage 8A ConditionGeneResponseRegistry.
        network_id: Optional network ID.
        network_name: Optional network title.
        min_conditions: Minimum required perturbation contrasts (default 8).
        min_chemicals: Minimum required distinct chemicals (default 4).
        method: Correlation method ('spearman' or 'pearson').
        metric_column: Differential statistic used for co-response ('log2_fold_change' or 'wald_statistic').
        min_abs_correlation: Minimum correlation cutoff to retain an edge (default 0.7).
        max_fdr: Maximum Benjamini-Hochberg FDR threshold (default 0.05).
        max_evaluated_genes: Optional limit on highest-variance genes to evaluate.
        training_chemicals: Optional set of compound IDs to enforce AI/ML training partition safety.
        partition_assignment: Safety partition flag.

    Returns:
        tuple of (BiologicalNetwork, CoresponseEvaluationSummary)
    """
    # 1. Eligibility Gating
    is_eligible, n_cond, n_chem, reason = validate_cross_condition_eligibility(
        registry=registry,
        min_conditions=min_conditions,
        min_chemicals=min_chemicals,
        training_chemicals=training_chemicals,
    )

    all_genes = registry.unique_canonical_genes
    net_id = network_id or f"coresponse_{registry.registry_id}"
    net_name = network_name or f"Cross-Condition Co-Response Network ({registry.registry_id})"

    if not is_eligible:
        qc_stat = CoresponseQCStatus.INSUFFICIENT_CONDITIONS if n_cond < min_conditions else CoresponseQCStatus.INSUFFICIENT_CHEMICALS
        summary = CoresponseEvaluationSummary(
            status=qc_stat,
            is_eligible=False,
            evaluated_genes_count=len(all_genes),
            evaluated_conditions_count=n_cond,
            evaluated_chemicals_count=n_chem,
            evaluated_samples_count=0,
            inferred_edges_count=0,
            reasons=(reason,),
            scientific_limitations=reason,
        )
        empty_net = BiologicalNetwork(
            network_id=net_id,
            network_name=f"{net_name} - INSUFFICIENT EVIDENCE",
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            source_database=registry.registry_id,
            database_version="empirical_v1",
            qc_status=qc_stat,
            scientific_limitations=reason,
            provenance={
                "registry_id": registry.registry_id,
                "conditions_found": n_cond,
                "chemicals_found": n_chem,
                "status": str(qc_stat),
            },
        )
        return empty_net, summary

    # 2. Build genes x contrasts matrix
    records = registry.records
    if training_chemicals is not None:
        records = tuple(r for r in records if r.compound_id in training_chemicals or r.compound_name in training_chemicals)

    rows = []
    gene_symbol_map: dict[str, str] = {}
    for r in records:
        val = getattr(r, metric_column, None)
        if val is None or np.isnan(val):
            continue
        rows.append({
            "gene_id": r.canonical_ensembl_gene_id,
            "contrast_id": r.contrast_id,
            "metric": float(val),
        })
        if r.gene_symbol:
            gene_symbol_map[r.canonical_ensembl_gene_id] = r.gene_symbol

    if not rows:
        reason_no_recs = f"No valid numerical entries found for metric '{metric_column}'."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            is_eligible=False,
            evaluated_genes_count=0,
            evaluated_conditions_count=n_cond,
            evaluated_chemicals_count=n_chem,
            evaluated_samples_count=0,
            inferred_edges_count=0,
            reasons=(reason_no_recs,),
            scientific_limitations=reason_no_recs,
        )
        empty_net = BiologicalNetwork(
            network_id=net_id,
            network_name=net_name,
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            qc_status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            scientific_limitations=reason_no_recs,
        )
        return empty_net, summary

    df_records = pd.DataFrame(rows)
    # Pivot: rows = gene_id, cols = contrast_id
    pivoted = df_records.pivot(index="gene_id", columns="contrast_id", values="metric").dropna()
    if pivoted.shape[1] < min_conditions:
        reason_drop = f"After dropping genes with missing values across contrasts, only {pivoted.shape[1]} contrast(s) remain."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.INSUFFICIENT_CONDITIONS,
            is_eligible=False,
            evaluated_genes_count=pivoted.shape[0],
            evaluated_conditions_count=pivoted.shape[1],
            evaluated_chemicals_count=n_chem,
            evaluated_samples_count=0,
            inferred_edges_count=0,
            reasons=(reason_drop,),
            scientific_limitations=reason_drop,
        )
        empty_net = BiologicalNetwork(
            network_id=net_id,
            network_name=net_name,
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            qc_status=CoresponseQCStatus.INSUFFICIENT_CONDITIONS,
            scientific_limitations=reason_drop,
        )
        return empty_net, summary

    # 3. Variance filtering
    variances = pivoted.var(axis=1)
    valid_genes = variances[variances > 1e-9].index.tolist()
    if not valid_genes:
        reason_var = "All evaluated genes exhibit zero response variance across perturbation contrasts."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            is_eligible=False,
            evaluated_genes_count=0,
            evaluated_conditions_count=n_cond,
            evaluated_chemicals_count=n_chem,
            evaluated_samples_count=0,
            inferred_edges_count=0,
            reasons=(reason_var,),
            scientific_limitations=reason_var,
        )
        empty_net = BiologicalNetwork(
            network_id=net_id,
            network_name=net_name,
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            qc_status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            scientific_limitations=reason_var,
        )
        return empty_net, summary

    if max_evaluated_genes and len(valid_genes) > max_evaluated_genes:
        top_genes = variances.loc[valid_genes].nlargest(max_evaluated_genes).index.tolist()
        mat_sub = pivoted.loc[top_genes]
    else:
        top_genes = valid_genes
        mat_sub = pivoted.loc[top_genes]

    mat = mat_sub.values  # Shape: (G, C)
    num_genes = mat.shape[0]
    n_contrasts = mat.shape[1]

    if method.lower() == "spearman":
        ranked_mat = np.apply_along_axis(stats.rankdata, 1, mat)
        corr_matrix = np.corrcoef(ranked_mat)
    else:
        corr_matrix = np.corrcoef(mat)

    # 4. Pairwise t-test and FDR
    df_deg = n_contrasts - 2
    raw_edges = []
    p_vals = []

    for i in range(num_genes):
        for j in range(i + 1, num_genes):
            r = corr_matrix[i, j]
            if np.isnan(r):
                continue
            abs_r = abs(r)
            if abs_r < min_abs_correlation:
                continue

            clipped_r2 = min(0.999999, r * r)
            t_stat = r * math.sqrt(df_deg / (1.0 - clipped_r2))
            p_val = float(stats.t.sf(abs(t_stat), df=df_deg) * 2.0)

            raw_edges.append((top_genes[i], top_genes[j], r, p_val))
            p_vals.append(p_val)

    if not raw_edges:
        reason_no_sig = f"No gene pairs satisfied minimum co-response correlation cutoff (|r| >= {min_abs_correlation})."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.VALID,
            is_eligible=True,
            evaluated_genes_count=num_genes,
            evaluated_conditions_count=n_contrasts,
            evaluated_chemicals_count=n_chem,
            evaluated_samples_count=0,
            inferred_edges_count=0,
            reasons=(reason_no_sig,),
            scientific_limitations=reason_no_sig,
        )
        empty_net = BiologicalNetwork(
            network_id=net_id,
            network_name=net_name,
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            qc_status=CoresponseQCStatus.VALID,
            scientific_limitations=reason_no_sig,
        )
        return empty_net, summary

    fdr_vals = compute_benjamini_hochberg_fdr(p_vals)
    final_edges: list[GeneNetworkEdge] = []
    distinct_compounds = tuple(sorted(list({r.compound_name for r in records})))
    distinct_studies = tuple(sorted(list({r.study_id for r in records})))

    for (g1, g2, r, p_val), fdr in zip(raw_edges, fdr_vals):
        if fdr > max_fdr:
            continue

        sym1 = gene_symbol_map.get(g1)
        sym2 = gene_symbol_map.get(g2)
        edge_id = derive_edge_id(g1, g2, NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE, registry.registry_id)

        edge = GeneNetworkEdge(
            edge_id=edge_id,
            gene_a_canonical_id=g1,
            gene_b_canonical_id=g2,
            gene_a_symbol=sym1,
            gene_b_symbol=sym2,
            relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
            interaction_nature=InteractionNature.STATISTICAL_CO_RESPONSE,
            evidence_source=registry.registry_id,
            evidence_version="empirical_v1",
            evidence_category=EvidenceCategory.EMPIRICAL_CORESPONSE,
            confidence_score=float(abs(r)),
            statistical_metric_name=f"{method.lower()}_coresponse",
            statistical_estimate=float(r),
            p_value=float(p_val),
            adjusted_p_value=float(fdr),
            condition_count=n_contrasts,
            study_context=distinct_studies,
            compound_context=distinct_compounds,
            partition_assignment=partition_assignment,
            qc_status=CoresponseQCStatus.VALID,
            is_direct_binding=False,
            is_chemical_specific=True,  # Empirical co-response across chemical perturbation conditions
            causality_claim="COORDINATED_TRANSCRIPTIONAL_RESPONSE",
            scientific_limitations=(
                f"Statistically significant cross-condition co-response ({method}, r={r:.3f}, FDR={fdr:.3e}) "
                f"across {n_contrasts} perturbation conditions. Does NOT imply direct physical binding or shared transcriptional regulation."
            ),
            provenance_notes=f"Registry: {registry.registry_id}; Contrasts: {n_contrasts}; Chemicals: {n_chem}",
        )
        final_edges.append(edge)

    summary = CoresponseEvaluationSummary(
        status=CoresponseQCStatus.VALID,
        is_eligible=True,
        evaluated_genes_count=num_genes,
        evaluated_conditions_count=n_contrasts,
        evaluated_chemicals_count=n_chem,
        evaluated_samples_count=0,
        inferred_edges_count=len(final_edges),
        reasons=(f"Successfully identified {len(final_edges)} cross-condition co-response edges (FDR <= {max_fdr}).",),
        scientific_limitations="Empirical cross-condition response coordination; non-causal observational association.",
    )

    net = BiologicalNetwork(
        network_id=net_id,
        network_name=net_name,
        relationship_type=NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE,
        edges=tuple(final_edges),
        species="Homo sapiens",
        taxon_id=9606,
        source_database=registry.registry_id,
        database_version="empirical_v1",
        qc_status=CoresponseQCStatus.VALID,
        scientific_limitations=summary.scientific_limitations,
        provenance={
            "registry_id": registry.registry_id,
            "method": method,
            "metric_column": metric_column,
            "min_abs_correlation": min_abs_correlation,
            "max_fdr": max_fdr,
            "condition_count": n_contrasts,
            "chemical_count": n_chem,
            "edge_count": len(final_edges),
        },
    )
    return net, summary
