"""End-to-end orchestration pipeline for Stage 8C Pathway and Signaling Analysis."""

from pathlib import Path
from typing import Sequence

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
)
from src.gene_harmonization.reference import GeneReference

from .mapping import EnsemblToEntrezAdapter
from .models import (
    PathwayAnalysisProvenance,
    PathwayAnalysisQC,
    PathwayEnrichmentDataset,
    PathwayEnrichmentResult,
    PathwayGeneMappingResult,
    RankedGeneSet,
    RetainedNonEntrezGeneFeature,
    StatisticallyEligibleGeneFeature,
)
from .provenance import build_pathway_analysis_provenance
from .qc import compute_pathway_analysis_qc
from .r_bridge import execute_preranked_gsea_r
from .ranking import GeneRankingBuilder


def run_pathway_enrichment_for_contrast(
    contrast: DifferentialExpressionContrastResult,
    gene_reference: GeneReference | None = None,
    r_binary_path: str = "Rscript",
    min_size: int = 15,
    max_size: int = 500,
    pvalue_cutoff: float = 1.0,
    seed: int = 42,
    gsea_param: float = 1.0,
    score_type: str = "std",
    nperm_simple: int = 1000,
    timeout_seconds: int = 600,
    deduplicate_entrez: bool = True,
) -> PathwayEnrichmentDataset:
    """Run full Stage 8C pathway and signaling analysis for a single differential expression contrast.

    Steps:
    1. Map canonical Ensembl IDs to Entrez Gene IDs using governed local reference.
    2. Build statistically qualified, deterministically ordered ranked gene set.
    3. Preserve the full statistically eligible molecular universe (32,140 genes) and retained non-Entrez genes.
    4. Execute preranked GSEA in R using clusterProfiler / fgsea for GO BP and KEGG.
    5. Compute exhaustive quality control accounting and metrics.
    6. Construct complete upstream-to-downstream provenance audit record.
    7. Return governed PathwayEnrichmentDataset.
    """
    mapping_adapter = EnsemblToEntrezAdapter(gene_reference=gene_reference)
    ranking_builder = GeneRankingBuilder(mapping_adapter=mapping_adapter)

    ranking_res = ranking_builder.build_ranking(
        contrast=contrast,
        deduplicate_entrez=deduplicate_entrez,
    )
    ranked_genes = ranking_res.ranked_gene_set
    mapping_result = ranking_res.mapping_result
    excluded_records = ranking_res.excluded_records
    stat_universe = ranking_res.statistically_eligible_universe
    retained_non_entrez = ranking_res.retained_non_entrez_genes

    # Execute GSEA via R subprocess bridge
    runner_result = execute_preranked_gsea_r(
        ranked_genes=ranked_genes,
        r_binary_path=r_binary_path,
        min_size=min_size,
        max_size=max_size,
        pvalue_cutoff=pvalue_cutoff,
        seed=seed,
        gsea_param=gsea_param,
        score_type=score_type,
        nperm_simple=nperm_simple,
        timeout_seconds=timeout_seconds,
    )

    if runner_result.status != "OK":
        raise RuntimeError(
            f"Stage 8C GSEA execution failed ({runner_result.status}): {runner_result.error_message}"
        )

    qc = compute_pathway_analysis_qc(
        stage7_gene_count=len(contrast.gene_results),
        mapping_result=mapping_result,
        ranked_set=ranked_genes,
        excluded_records=excluded_records,
        go_results=runner_result.go_results,
        kegg_results=runner_result.kegg_results,
        contrast=contrast,
        statistically_eligible_universe=stat_universe,
        retained_non_entrez_genes=retained_non_entrez,
    )

    provenance = build_pathway_analysis_provenance(
        contrast=contrast,
        ranking_metric=ranked_genes.ranking_metric_name,
        gene_reference_version=mapping_result.reference_version,
        pathway_databases=("GO_BP", "KEGG"),
        software_packages=runner_result.packages,
        gsea_parameters=runner_result.gsea_parameters,
    )

    return PathwayEnrichmentDataset(
        contrast_id=contrast.contrast_id,
        study_id=contrast.study_id,
        go_results=runner_result.go_results,
        kegg_results=runner_result.kegg_results,
        ranked_genes=ranked_genes,
        statistically_eligible_universe=stat_universe,
        retained_non_entrez_genes=retained_non_entrez,
        qc=qc,
        provenance=provenance,
        mapping_result=mapping_result,
    )


def export_pathway_results_to_csv(
    results: Sequence[PathwayEnrichmentResult],
    output_path: Path | str,
) -> None:
    """Export pathway enrichment results to a transparent tabular CSV file."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "contrast_id",
        "source_database",
        "pathway_id",
        "pathway_name",
        "enrichment_score",
        "normalized_enrichment_score",
        "p_value",
        "adjusted_p_value",
        "pathway_size",
        "direction",
        "is_significant",
        "leading_edge_symbols",
        "leading_edge_gene_ids",
        "ranking_metric",
        "analysis_version",
    ]

    lines = [",".join(header)]
    for r in results:
        symbols_str = ";".join(r.leading_edge_symbols)
        ids_str = ";".join(r.leading_edge_gene_ids)
        name_clean = f'"{r.pathway_name}"' if "," in r.pathway_name else r.pathway_name
        row = [
            r.contrast_id,
            r.source_database,
            r.pathway_id,
            name_clean,
            f"{r.enrichment_score:.6f}",
            f"{r.normalized_enrichment_score:.6f}",
            f"{r.p_value:.6e}",
            f"{r.adjusted_p_value:.6e}",
            str(r.pathway_size),
            r.direction,
            str(r.is_significant),
            f'"{symbols_str}"' if "," in symbols_str or ";" in symbols_str else symbols_str,
            f'"{ids_str}"' if "," in ids_str or ";" in ids_str else ids_str,
            r.ranking_metric,
            r.analysis_version,
        ]
        lines.append(",".join(row))

    with open(p, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def export_retained_non_entrez_genes_to_csv(
    genes: Sequence[RetainedNonEntrezGeneFeature],
    output_path: Path | str,
) -> None:
    """Export retained non-Entrez statistically eligible genes to a tabular CSV file."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "canonical_gene_id",
        "approved_symbol",
        "gene_biotype",
        "wald_statistic",
        "log2_fold_change",
        "lfc_standard_error",
        "p_value",
        "adjusted_p_value_bh",
        "reason_not_in_kegg",
        "go_annotation_status",
        "future_feature_eligibility",
    ]

    lines = [",".join(header)]
    for g in genes:
        row = [
            g.canonical_gene_id,
            g.approved_symbol or "",
            g.gene_biotype or "",
            f"{g.wald_statistic:.6f}" if g.wald_statistic is not None else "",
            f"{g.log2_fold_change:.6f}" if g.log2_fold_change is not None else "",
            f"{g.lfc_standard_error:.6f}" if g.lfc_standard_error is not None else "",
            f"{g.p_value:.6e}" if g.p_value is not None else "",
            f"{g.adjusted_p_value_bh:.6e}" if g.adjusted_p_value_bh is not None else "",
            g.reason_not_in_kegg,
            g.go_annotation_status,
            str(g.future_feature_eligibility),
        ]
        lines.append(",".join(row))

    with open(p, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def export_statistically_eligible_universe_to_csv(
    universe: Sequence[StatisticallyEligibleGeneFeature],
    output_path: Path | str,
) -> None:
    """Export the complete preserved statistically eligible molecular universe to CSV."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "canonical_gene_id",
        "approved_symbol",
        "gene_biotype",
        "entrez_gene_id",
        "base_mean",
        "log2_fold_change",
        "lfc_standard_error",
        "wald_statistic",
        "p_value",
        "adjusted_p_value_bh",
        "result_status",
        "annotation_status",
        "statistical_eligibility",
    ]

    lines = [",".join(header)]
    for g in universe:
        row = [
            g.canonical_gene_id,
            g.approved_symbol or "",
            g.gene_biotype or "",
            g.entrez_gene_id or "",
            f"{g.base_mean:.6f}",
            f"{g.log2_fold_change:.6f}" if g.log2_fold_change is not None else "",
            f"{g.lfc_standard_error:.6f}" if g.lfc_standard_error is not None else "",
            f"{g.wald_statistic:.6f}" if g.wald_statistic is not None else "",
            f"{g.p_value:.6e}" if g.p_value is not None else "",
            f"{g.adjusted_p_value_bh:.6e}" if g.adjusted_p_value_bh is not None else "",
            g.result_status,
            g.annotation_status,
            str(g.statistical_eligibility),
        ]
        lines.append(",".join(row))

    with open(p, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
