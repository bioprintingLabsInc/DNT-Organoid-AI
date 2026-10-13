"""Quality control computation and integrity verification for Stage 8C."""

from typing import Sequence

from .models import (
    PathwayAnalysisQC,
    PathwayEnrichmentResult,
    PathwayGeneMappingResult,
    RankedGeneSet,
)
from .ranking import ExcludedRankingGene


def compute_pathway_analysis_qc(
    stage7_gene_count: int,
    mapping_result: PathwayGeneMappingResult,
    ranked_set: RankedGeneSet,
    excluded_records: Sequence[ExcludedRankingGene],
    go_results: Sequence[PathwayEnrichmentResult],
    kegg_results: Sequence[PathwayEnrichmentResult],
) -> PathwayAnalysisQC:
    """Compute exhaustive Stage 8C pathway analysis quality control report."""
    total = stage7_gene_count if stage7_gene_count > 0 else 1

    valid_stat_count = total - mapping_result.invalid_stat_count
    valid_stat_pct = round((valid_stat_count / total) * 100.0, 2)
    invalid_stat_pct = round((mapping_result.invalid_stat_count / total) * 100.0, 2)

    mapped_pct = round((mapping_result.mapped_count / total) * 100.0, 2)
    unmapped_pct = round((mapping_result.unmapped_count / total) * 100.0, 2)

    genes_go = len(ranked_set.ranked_genes)
    genes_kegg = len(ranked_set.ranked_genes)

    go_tested = len(go_results)
    kegg_tested = len(kegg_results)

    go_sig = sum(1 for r in go_results if r.is_significant)
    go_sig_pct = round((go_sig / go_tested) * 100.0, 2) if go_tested > 0 else 0.0

    kegg_sig = sum(1 for r in kegg_results if r.is_significant)
    kegg_sig_pct = round((kegg_sig / kegg_tested) * 100.0, 2) if kegg_tested > 0 else 0.0

    notes: list[str] = []
    if genes_go == 0:
        qc_status = "FAIL"
        notes.append("Zero genes eligible for GSEA analysis.")
    elif go_tested == 0 and kegg_tested == 0:
        qc_status = "FAIL"
        notes.append("Zero gene sets evaluated in GSEA.")
    elif unmapped_pct > 50.0:
        qc_status = "WARNING"
        notes.append(f"Unmapped gene rate is {unmapped_pct}% (expected due to non-coding Ensembl genes lacking Entrez IDs).")
    else:
        qc_status = "PASS"
        notes.append("GSEA input ranking and pathway evaluation completed successfully.")

    if mapping_result.duplicate_entrez_count > 0:
        notes.append(
            f"Resolved {mapping_result.duplicate_entrez_count} duplicate Entrez mappings using {ranked_set.duplicate_resolution_rule}."
        )

    return PathwayAnalysisQC(
        stage7_evaluated_genes=stage7_gene_count,
        valid_ranking_statistic_genes=valid_stat_count,
        valid_ranking_statistic_pct=valid_stat_pct,
        invalid_statistical_result_genes=mapping_result.invalid_stat_count,
        invalid_statistical_result_pct=invalid_stat_pct,
        mapped_entrez_genes=mapping_result.mapped_count,
        mapped_entrez_pct=mapped_pct,
        unmapped_ensembl_genes=mapping_result.unmapped_count,
        unmapped_ensembl_pct=unmapped_pct,
        one_to_many_mappings=mapping_result.one_to_many_count,
        duplicate_entrez_mappings=mapping_result.duplicate_entrez_count,
        genes_entering_go_gsea=genes_go,
        genes_entering_kegg_gsea=genes_kegg,
        go_gene_sets_tested=go_tested,
        kegg_gene_sets_tested=kegg_tested,
        go_significant_pathways_count=go_sig,
        go_significant_pathways_pct=go_sig_pct,
        kegg_significant_pathways_count=kegg_sig,
        kegg_significant_pathways_pct=kegg_sig_pct,
        qc_status=qc_status,
        qc_notes=tuple(notes),
    )
