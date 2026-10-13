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
    timeout_seconds: int = 600,
    deduplicate_entrez: bool = True,
) -> PathwayEnrichmentDataset:
    """Run full Stage 8C pathway and signaling analysis for a single differential expression contrast.

    Steps:
    1. Map canonical Ensembl IDs to Entrez Gene IDs using governed local reference.
    2. Build statistically qualified, deterministically ordered ranked gene set.
    3. Execute preranked GSEA in R using clusterProfiler / fgsea for GO BP and KEGG.
    4. Compute exhaustive quality control accounting and metrics.
    5. Construct complete upstream-to-downstream provenance audit record.
    6. Return governed PathwayEnrichmentDataset.
    """
    mapping_adapter = EnsemblToEntrezAdapter(gene_reference=gene_reference)
    ranking_builder = GeneRankingBuilder(mapping_adapter=mapping_adapter)

    ranked_genes, mapping_result, excluded_records = ranking_builder.build_ranking(
        contrast=contrast,
        deduplicate_entrez=deduplicate_entrez,
    )

    # Execute GSEA via R subprocess bridge
    runner_result = execute_preranked_gsea_r(
        ranked_genes=ranked_genes,
        r_binary_path=r_binary_path,
        min_size=min_size,
        max_size=max_size,
        pvalue_cutoff=pvalue_cutoff,
        seed=seed,
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
    )

    provenance = build_pathway_analysis_provenance(
        contrast=contrast,
        ranking_metric=ranked_genes.ranking_metric_name,
        gene_reference_version=mapping_result.reference_version,
        pathway_databases=("GO_BP", "KEGG"),
        software_packages=runner_result.packages,
    )

    return PathwayEnrichmentDataset(
        contrast_id=contrast.contrast_id,
        study_id=contrast.study_id,
        go_results=runner_result.go_results,
        kegg_results=runner_result.kegg_results,
        ranked_genes=ranked_genes,
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

    p.write_text("\n".join(lines), encoding="utf-8")
