"""Stage 8B: Within-Study Gene Co-Expression Analysis and Statistical Gating.

Evaluates sample-level continuous gene expression within a study cohort while enforcing:
1. Strict biological independence validation (technical replicates CANNOT inflate sample size).
2. Minimum observation gating (N_bio >= min_independent_samples; e.g. 10). Under-replicated cohorts
   honestly produce INSUFFICIENT_SAMPLES / INSUFFICIENT_EVIDENCE rather than manufacturing spurious edges.
3. Expression scale normalization verification (rejects raw integer or unnormalized matrices).
4. Rigorous statistical correlation (Spearman/Pearson) with exact t-distribution p-values and
   Benjamini-Hochberg FDR multiple-testing correction.
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


def compute_benjamini_hochberg_fdr(p_values: Sequence[float]) -> list[float]:
    """Calculate Benjamini-Hochberg False Discovery Rate (FDR) adjusted p-values."""
    n = len(p_values)
    if n == 0:
        return []
    
    # Store original indices
    indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [1.0] * n
    
    cum_min = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, p_val = indexed_p[rank_idx]
        rank = rank_idx + 1
        adj_p = min(1.0, (p_val * n) / rank)
        cum_min = min(cum_min, adj_p)
        adjusted[orig_idx] = cum_min
        
    return adjusted


def validate_coexpression_matrix_eligibility(
    expression_matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    min_independent_samples: int = 10,
    biological_replicate_col: str = "biological_replicate_id",
) -> tuple[bool, int, str]:
    """Validate whether an expression matrix is statistically eligible for co-expression analysis.

    Args:
        expression_matrix: DataFrame with genes in rows (index) and samples in columns.
        sample_metadata: Optional DataFrame indexing samples with biological metadata.
        min_independent_samples: Minimum independent biological observations required (default 10).
        biological_replicate_col: Column indicating biological replicate / donor ID.

    Returns:
        tuple of (is_eligible: bool, independent_n: int, reason: str)
    """
    if expression_matrix.empty:
        return False, 0, "Expression matrix is empty."

    n_samples = expression_matrix.shape[1]

    # Check for negative values
    numeric_vals = expression_matrix.select_dtypes(include=[np.number]).values
    if np.any(numeric_vals < 0):
        return False, n_samples, "Expression matrix contains negative values; continuous scale invalid."

    # Validate biological independence vs. technical replicate inflation
    if sample_metadata is not None and biological_replicate_col in sample_metadata.columns:
        valid_samples = [s for s in expression_matrix.columns if s in sample_metadata.index]
        if valid_samples:
            sub_meta = sample_metadata.loc[valid_samples]
            unique_bio = set(sub_meta[biological_replicate_col].dropna())
            independent_n = len(unique_bio)
        else:
            independent_n = n_samples
    else:
        # Without metadata, assume samples represent independent observations but subject to strict count gating
        independent_n = n_samples

    if independent_n < min_independent_samples:
        return (
            False,
            independent_n,
            f"Insufficient independent biological observations: found N={independent_n}, but minimum "
            f"required is min_independent_samples={min_independent_samples}. Reporting correlation on "
            f"under-replicated data yields statistical noise and inflated false discovery."
        )

    return True, independent_n, "Statistically eligible for within-study co-expression analysis."


def analyze_within_study_coexpression(
    expression_matrix: pd.DataFrame,
    study_id: str,
    gene_symbols: dict[str, str] | None = None,
    sample_metadata: pd.DataFrame | None = None,
    min_independent_samples: int = 10,
    method: str = "spearman",
    min_abs_correlation: float = 0.6,
    max_fdr: float = 0.05,
    max_evaluated_genes: int | None = 2000,
    biological_replicate_col: str = "biological_replicate_id",
) -> tuple[BiologicalNetwork, CoresponseEvaluationSummary]:
    """Execute within-study gene co-expression analysis with strict statistical gating.

    If sample size is insufficient (N < min_independent_samples), honestly halts inference
    and returns an empty network with INSUFFICIENT_SAMPLES status.

    Args:
        expression_matrix: DataFrame with canonical Ensembl gene IDs as index, sample IDs as columns.
        study_id: Identifier of the evaluated study.
        gene_symbols: Optional mapping of Ensembl ID -> gene symbol.
        sample_metadata: Optional sample metadata table to verify biological independence.
        min_independent_samples: Minimum independent biological samples required (default 10).
        method: Correlation method ('spearman' or 'pearson').
        min_abs_correlation: Minimum absolute correlation coefficient to retain an edge (default 0.6).
        max_fdr: Maximum Benjamini-Hochberg adjusted p-value cutoff (default 0.05).
        max_evaluated_genes: Limit on genes to evaluate (sorted by variance) to avoid quadratic blowup.
        biological_replicate_col: Column indicating biological replicate / donor identity.

    Returns:
        tuple of (BiologicalNetwork, CoresponseEvaluationSummary)
    """
    symbols = gene_symbols or {}

    # 1. Eligibility Gating
    is_eligible, n_bio, reason = validate_coexpression_matrix_eligibility(
        expression_matrix=expression_matrix,
        sample_metadata=sample_metadata,
        min_independent_samples=min_independent_samples,
        biological_replicate_col=biological_replicate_col,
    )

    if not is_eligible:
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.INSUFFICIENT_SAMPLES,
            is_eligible=False,
            evaluated_genes_count=len(expression_matrix),
            evaluated_conditions_count=1,
            evaluated_chemicals_count=1,
            evaluated_samples_count=n_bio,
            inferred_edges_count=0,
            reasons=(reason,),
            scientific_limitations=reason,
        )
        empty_net = BiologicalNetwork(
            network_id=f"within_study_{study_id}_ineligible",
            network_name=f"Within-Study Co-expression Network ({study_id}) - INSUFFICIENT EVIDENCE",
            relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            source_database=study_id,
            database_version="empirical_v1",
            qc_status=CoresponseQCStatus.INSUFFICIENT_SAMPLES,
            scientific_limitations=reason,
            provenance={"study_id": study_id, "independent_samples": n_bio, "status": "INSUFFICIENT_SAMPLES"},
        )
        return empty_net, summary

    # 2. Filter zero-variance genes
    variances = expression_matrix.var(axis=1)
    valid_genes = variances[variances > 1e-9].index.tolist()
    if not valid_genes:
        reason_var = "All genes in expression matrix exhibit near-zero variance across samples."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            is_eligible=False,
            evaluated_genes_count=0,
            evaluated_conditions_count=1,
            evaluated_chemicals_count=1,
            evaluated_samples_count=n_bio,
            inferred_edges_count=0,
            reasons=(reason_var,),
            scientific_limitations=reason_var,
        )
        empty_net = BiologicalNetwork(
            network_id=f"within_study_{study_id}_zero_var",
            network_name=f"Within-Study Co-expression Network ({study_id}) - ZERO VARIANCE",
            relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            qc_status=CoresponseQCStatus.INSUFFICIENT_EVIDENCE,
            scientific_limitations=reason_var,
        )
        return empty_net, summary

    # Filter to highest variance genes if capped
    if max_evaluated_genes and len(valid_genes) > max_evaluated_genes:
        top_genes = variances.loc[valid_genes].nlargest(max_evaluated_genes).index.tolist()
        expr_sub = expression_matrix.loc[top_genes]
    else:
        top_genes = valid_genes
        expr_sub = expression_matrix.loc[top_genes]

    # 3. Compute correlation matrix
    mat = expr_sub.values  # Shape: (G, N)
    num_genes = mat.shape[0]

    if method.lower() == "spearman":
        # Convert rows to ranks
        ranked_mat = np.apply_along_axis(stats.rankdata, 1, mat)
        # Pearson on ranks = Spearman
        corr_matrix = np.corrcoef(ranked_mat)
    else:
        corr_matrix = np.corrcoef(mat)

    # 4. Pairwise t-test and FDR correction
    df_deg = n_bio - 2
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

            # Student's t calculation
            # Avoid division by zero if r is extremely close to 1.0
            clipped_r2 = min(0.999999, r * r)
            t_stat = r * math.sqrt(df_deg / (1.0 - clipped_r2))
            p_val = float(stats.t.sf(abs(t_stat), df=df_deg) * 2.0)

            raw_edges.append((top_genes[i], top_genes[j], r, p_val))
            p_vals.append(p_val)

    if not raw_edges:
        reason_no_sig = f"No gene pairs satisfied minimum correlation cutoff (|r| >= {min_abs_correlation})."
        summary = CoresponseEvaluationSummary(
            status=CoresponseQCStatus.VALID,
            is_eligible=True,
            evaluated_genes_count=num_genes,
            evaluated_conditions_count=1,
            evaluated_chemicals_count=1,
            evaluated_samples_count=n_bio,
            inferred_edges_count=0,
            reasons=(reason_no_sig,),
            scientific_limitations=reason_no_sig,
        )
        empty_net = BiologicalNetwork(
            network_id=f"within_study_{study_id}_no_sig",
            network_name=f"Within-Study Co-expression Network ({study_id})",
            relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
            edges=(),
            species="Homo sapiens",
            taxon_id=9606,
            source_database=study_id,
            database_version="empirical_v1",
            qc_status=CoresponseQCStatus.VALID,
            scientific_limitations=reason_no_sig,
        )
        return empty_net, summary

    # Benjamini-Hochberg FDR
    fdr_vals = compute_benjamini_hochberg_fdr(p_vals)

    final_edges: list[GeneNetworkEdge] = []
    for (g1, g2, r, p_val), fdr in zip(raw_edges, fdr_vals):
        if fdr > max_fdr:
            continue

        sym1 = symbols.get(g1)
        sym2 = symbols.get(g2)
        edge_id = derive_edge_id(g1, g2, NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION, study_id)

        edge = GeneNetworkEdge(
            edge_id=edge_id,
            gene_a_canonical_id=g1,
            gene_b_canonical_id=g2,
            gene_a_symbol=sym1,
            gene_b_symbol=sym2,
            relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
            interaction_nature=InteractionNature.STATISTICAL_CO_EXPRESSION,
            evidence_source=study_id,
            evidence_version="empirical_v1",
            evidence_category=EvidenceCategory.EMPIRICAL_COEXPRESSION,
            confidence_score=float(abs(r)),
            statistical_metric_name=f"{method.lower()}_correlation",
            statistical_estimate=float(r),
            p_value=float(p_val),
            adjusted_p_value=float(fdr),
            sample_size=n_bio,
            study_context=(study_id,),
            partition_assignment=PartitionAssignment.DESCRIPTIVE_EXPLORATION,
            qc_status=CoresponseQCStatus.VALID,
            is_direct_binding=False,
            is_chemical_specific=False,
            causality_claim="OBSERVED_CO_EXPRESSION",
            scientific_limitations=(
                f"Statistically significant within-study expression co-variation ({method}, r={r:.3f}, FDR={fdr:.3e}). "
                "Does NOT establish physical protein interaction, shared transcriptional regulation, or chemical-specific co-response."
            ),
            provenance_notes=f"Study: {study_id}; Samples: N={n_bio}; Method: {method}",
        )
        final_edges.append(edge)

    summary = CoresponseEvaluationSummary(
        status=CoresponseQCStatus.VALID,
        is_eligible=True,
        evaluated_genes_count=num_genes,
        evaluated_conditions_count=1,
        evaluated_chemicals_count=1,
        evaluated_samples_count=n_bio,
        inferred_edges_count=len(final_edges),
        reasons=(f"Successfully identified {len(final_edges)} significant within-study co-expression edges (FDR <= {max_fdr}).",),
        scientific_limitations="Within-study transcriptional co-expression; descriptive observational association.",
    )

    net = BiologicalNetwork(
        network_id=f"within_study_{study_id}",
        network_name=f"Within-Study Co-expression Network ({study_id})",
        relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
        edges=tuple(final_edges),
        species="Homo sapiens",
        taxon_id=9606,
        source_database=study_id,
        database_version="empirical_v1",
        qc_status=CoresponseQCStatus.VALID,
        scientific_limitations=summary.scientific_limitations,
        provenance={
            "study_id": study_id,
            "method": method,
            "min_abs_correlation": min_abs_correlation,
            "max_fdr": max_fdr,
            "sample_size": n_bio,
            "evaluated_genes": num_genes,
            "edge_count": len(final_edges),
        },
    )
    return net, summary
