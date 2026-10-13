"""Integration tests for Stage 8C headless R GSEA runner and full contrast pipeline."""

import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.pathway.pipeline import run_pathway_enrichment_for_contrast
from src.features.pathway.r_bridge import execute_preranked_gsea_r
from src.features.pathway.ranking import GeneRankingBuilder


@pytest.fixture
def mini_ranked_contrast() -> DifferentialExpressionContrastResult:
    """Fixture providing a synthetic contrast with 500 real canonical genes from Ensembl 112."""
    from src.features.pathway.mapping import EnsemblToEntrezAdapter

    adapter = EnsemblToEntrezAdapter()
    # Select first 400 protein coding genes with Entrez IDs
    gene_recs = [
        r for r in adapter.gene_reference.records
        if r.entrez_gene_id and "|" not in r.entrez_gene_id and r.gene_biotype == "protein_coding"
    ][:400]

    gene_results = []
    for idx, r in enumerate(gene_recs):
        # Create a spread of positive and negative wald statistics
        wald = 6.0 - (idx * 0.03)  # from +6.0 down to -6.0
        p_val = 1e-4 if abs(wald) > 2.0 else 0.2
        padj = 1e-3 if abs(wald) > 2.0 else 0.4
        gene_results.append(
            GeneDifferentialExpressionResult(
                canonical_gene_id=r.ensembl_gene_id,
                approved_symbol=r.gene_symbol,
                base_mean=100.0,
                log2_fold_change=wald * 0.5,
                lfc_standard_error=0.5,
                wald_statistic=wald,
                p_value=p_val,
                adjusted_p_value_bh=padj,
                result_status="OK",
            )
        )

    return DifferentialExpressionContrastResult(
        contrast_id="test_gsea_contrast_001",
        study_id="GSE_TEST",
        experiment_id="EXP_TEST",
        cohort_id="COHORT_TEST",
        relationship_id="REL_TEST",
        treatment_condition_id="TRT",
        matched_control_condition_ids=("CTL",),
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=("R1", "R2"),
        control_biological_replicate_ids=("R3", "R4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name="test_compound",
        agent_identifier="CID:9999",
        concentration_or_dose=50.0,
        concentration_or_dose_unit="uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 6",
        design_formula="~ condition",
        residual_degrees_of_freedom=2,
        size_factors_used={"S1": 1.0, "S2": 1.0, "S3": 1.0, "S4": 1.0},
        total_input_genes_count=len(gene_results),
        gene_results=tuple(gene_results),
    )


def test_gsea_runner_executes_and_returns_valid_schema(mini_ranked_contrast: DifferentialExpressionContrastResult) -> None:
    builder = GeneRankingBuilder()
    ranked_genes, mapping_res, _ = builder.build_ranking(mini_ranked_contrast)

    res = execute_preranked_gsea_r(
        ranked_genes=ranked_genes,
        min_size=10,
        max_size=500,
        pvalue_cutoff=1.0,
        seed=42,
        timeout_seconds=120,
    )

    assert res.status == "OK", f"GSEA runner failed: {res.error_message}"
    assert "clusterProfiler" in res.packages
    assert "org_Hs_eg_db" in res.packages
    assert "fgsea" in res.packages

    # Check GO results
    assert len(res.go_results) > 0, "Expected GO results to be returned"
    first_go = res.go_results[0]
    assert first_go.source_database == "GO_BP"
    assert first_go.pathway_id.startswith("GO:")
    assert len(first_go.pathway_name) > 0
    assert first_go.direction in ("POSITIVE", "NEGATIVE")
    assert 0.0 <= first_go.p_value <= 1.0
    assert 0.0 <= first_go.adjusted_p_value <= 1.0
    assert first_go.ranking_metric == "wald_statistic"

    # Check KEGG results
    assert len(res.kegg_results) > 0, "Expected KEGG results to be returned"
    first_kegg = res.kegg_results[0]
    assert first_kegg.source_database == "KEGG"
    assert first_kegg.pathway_id.startswith("hsa")
    assert first_kegg.direction in ("POSITIVE", "NEGATIVE")


def test_pipeline_end_to_end_orchestration(mini_ranked_contrast: DifferentialExpressionContrastResult) -> None:
    dataset = run_pathway_enrichment_for_contrast(
        contrast=mini_ranked_contrast,
        min_size=10,
        max_size=500,
        pvalue_cutoff=1.0,
        seed=42,
        timeout_seconds=180,
    )

    assert dataset.contrast_id == mini_ranked_contrast.contrast_id
    assert dataset.study_id == mini_ranked_contrast.study_id
    assert len(dataset.go_results) > 0
    assert len(dataset.kegg_results) > 0
    assert dataset.qc.qc_status in ("PASS", "WARNING")
    assert dataset.qc.stage7_evaluated_genes == len(mini_ranked_contrast.gene_results)
    assert dataset.provenance.agent_name == "test_compound"
    assert dataset.provenance.concentration == 50.0
