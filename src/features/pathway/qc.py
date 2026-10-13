"""Quality control computation and integrity verification for Stage 8C."""

import math
from typing import Sequence

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)

from .models import (
    PathwayAnalysisQC,
    PathwayEnrichmentResult,
    PathwayGeneMappingResult,
    RankedGeneSet,
    RetainedNonEntrezGeneFeature,
    StatisticallyEligibleGeneFeature,
)
from .ranking import ExcludedRankingGene


def compute_pathway_analysis_qc(
    stage7_gene_count: int,
    mapping_result: PathwayGeneMappingResult,
    ranked_set: RankedGeneSet,
    excluded_records: Sequence[ExcludedRankingGene],
    go_results: Sequence[PathwayEnrichmentResult],
    kegg_results: Sequence[PathwayEnrichmentResult],
    contrast: DifferentialExpressionContrastResult | None = None,
    statistically_eligible_universe: Sequence[StatisticallyEligibleGeneFeature] | None = None,
    retained_non_entrez_genes: Sequence[RetainedNonEntrezGeneFeature] | None = None,
) -> PathwayAnalysisQC:
    """Compute exhaustive Stage 8C pathway analysis quality control report."""
    total = stage7_gene_count if stage7_gene_count > 0 else 1

    if contrast is not None and contrast.gene_results:
        gene_results = contrast.gene_results
        total = len(gene_results)

        zero_count_genes = sum(
            1
            for g in gene_results
            if (
                g.base_mean == 0.0
                or (g.base_mean is not None and math.isclose(g.base_mean, 0.0, abs_tol=1e-9))
                or g.result_status == "ZERO_COUNTS"
            )
        )
        genes_with_valid_wald = sum(
            1
            for g in gene_results
            if g.wald_statistic is not None
            and not math.isnan(g.wald_statistic)
            and not math.isinf(g.wald_statistic)
        )
        cooks_outlier_genes = sum(
            1
            for g in gene_results
            if (
                g.wald_statistic is not None
                and not math.isnan(g.wald_statistic)
                and not math.isinf(g.wald_statistic)
                and g.result_status != "OK"
            )
        )
        statistically_eligible_genes = genes_with_valid_wald - cooks_outlier_genes
    elif statistically_eligible_universe is not None:
        statistically_eligible_genes = len(statistically_eligible_universe)
        outlier_count = sum(1 for e in excluded_records if e.result_status != "OK" and "Non-OK" in e.exclusion_reason)
        cooks_outlier_genes = outlier_count
        genes_with_valid_wald = statistically_eligible_genes + cooks_outlier_genes
        zero_count_genes = max(0, total - genes_with_valid_wald)
    else:
        # Fallback from excluded records and ranked set
        outlier_count = sum(1 for e in excluded_records if e.result_status != "OK")
        cooks_outlier_genes = outlier_count
        invalid_stat_count = sum(1 for e in excluded_records if "Invalid wald_statistic" in e.exclusion_reason)
        statistically_eligible_genes = total - invalid_stat_count - cooks_outlier_genes
        genes_with_valid_wald = statistically_eligible_genes + cooks_outlier_genes
        zero_count_genes = max(0, total - genes_with_valid_wald)

    if statistically_eligible_universe is not None:
        entrez_mapped_eligible_genes = sum(1 for g in statistically_eligible_universe if g.entrez_gene_id)
        entrez_unmapped_eligible_genes = sum(1 for g in statistically_eligible_universe if not g.entrez_gene_id)
    elif retained_non_entrez_genes is not None:
        entrez_unmapped_eligible_genes = len(retained_non_entrez_genes)
        entrez_mapped_eligible_genes = statistically_eligible_genes - entrez_unmapped_eligible_genes
    else:
        # Deduce from mapping_result or ranked_set
        dup_count = sum(1 for e in excluded_records if "Duplicate Entrez" in e.exclusion_reason)
        entrez_mapped_eligible_genes = len(ranked_set.ranked_genes) + dup_count
        entrez_unmapped_eligible_genes = max(0, statistically_eligible_genes - entrez_mapped_eligible_genes)

    dup_removed = sum(1 for e in excluded_records if "Duplicate Entrez" in e.exclusion_reason)
    if dup_removed == 0 and entrez_mapped_eligible_genes > len(ranked_set.ranked_genes):
        dup_removed = entrez_mapped_eligible_genes - len(ranked_set.ranked_genes)

    final_kegg_ranked = len(ranked_set.ranked_genes)
    final_go_ranked = len(ranked_set.ranked_genes)

    go_tested = len(go_results)
    kegg_tested = len(kegg_results)

    go_sig = sum(1 for r in go_results if r.is_significant)
    go_sig_pct = round((go_sig / go_tested) * 100.0, 2) if go_tested > 0 else 0.0

    kegg_sig = sum(1 for r in kegg_results if r.is_significant)
    kegg_sig_pct = round((kegg_sig / kegg_tested) * 100.0, 2) if kegg_tested > 0 else 0.0

    notes: list[str] = []
    if final_go_ranked == 0:
        qc_status = "FAIL"
        notes.append("Zero genes eligible for GSEA analysis.")
    elif go_tested == 0 and kegg_tested == 0:
        qc_status = "FAIL"
        notes.append("Zero gene sets evaluated in GSEA.")
    else:
        qc_status = "PASS"
        notes.append("GSEA input ranking and pathway evaluation completed successfully.")

    notes.append(
        f"Preserved STATISTICALLY_ELIGIBLE_MOLECULAR_UNIVERSE: {statistically_eligible_genes:,} genes "
        f"({round(statistically_eligible_genes / total * 100.0, 2)}% of evaluated canonical genes)."
    )
    if entrez_unmapped_eligible_genes > 0:
        notes.append(
            f"Retained {entrez_unmapped_eligible_genes:,} statistically valid non-Entrez genes "
            f"for cross-study molecular feature construction."
        )

    if dup_removed > 0:
        notes.append(
            f"Resolved {dup_removed} duplicate Entrez mappings using {ranked_set.duplicate_resolution_rule}."
        )

    return PathwayAnalysisQC(
        stage7_evaluated_genes=total,
        zero_count_genes=zero_count_genes,
        genes_with_valid_wald=genes_with_valid_wald,
        cooks_outlier_genes=cooks_outlier_genes,
        statistically_eligible_genes=statistically_eligible_genes,
        entrez_mapped_eligible_genes=entrez_mapped_eligible_genes,
        entrez_unmapped_eligible_genes=entrez_unmapped_eligible_genes,
        duplicate_entrez_entries_removed=dup_removed,
        final_kegg_ranked_genes=final_kegg_ranked,
        final_go_ranked_genes=final_go_ranked,
        go_gene_sets_tested=go_tested,
        kegg_gene_sets_tested=kegg_tested,
        go_significant_pathways_count=go_sig,
        go_significant_pathways_pct=go_sig_pct,
        kegg_significant_pathways_count=kegg_sig,
        kegg_significant_pathways_pct=kegg_sig_pct,
        qc_status=qc_status,
        qc_notes=tuple(notes),
    )
