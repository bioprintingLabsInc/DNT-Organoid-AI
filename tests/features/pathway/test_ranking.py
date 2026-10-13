"""Unit tests for Stage 8C GeneRankingBuilder."""

import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.pathway.mapping import EnsemblToEntrezAdapter
from src.features.pathway.ranking import GeneRankingBuilder
from src.gene_harmonization.reference import GeneReference, ReferenceRecord


@pytest.fixture
def test_contrast_and_ref() -> tuple[DifferentialExpressionContrastResult, GeneReference]:
    records = (
        ReferenceRecord("ENSG01", "GENE1", "101", "HGNC:1", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG02", "GENE2", "102", "HGNC:2", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG03", "GENE3", "103", "HGNC:3", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG04", "GENE4", "104", "HGNC:4", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG05", "GENE5", "105", "HGNC:5", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG06", "GENE6_DUP", "101", "HGNC:6", "protein_coding", (), "1", ""),  # Dup of 101
    )
    ref = GeneReference("r1", "v1", "GRCh38", "Ensembl", records, "chk1")

    gene_results = (
        # Gene 1: Positive wald, non-significant p-value (p=0.4) -> MUST BE INCLUDED
        GeneDifferentialExpressionResult("ENSG01", "GENE1", 10.0, 0.5, 0.5, 1.0, 0.4, 0.6, "OK"),
        # Gene 2: Strongly positive wald (wald=5.2)
        GeneDifferentialExpressionResult("ENSG02", "GENE2", 20.0, 2.6, 0.5, 5.2, 1e-5, 1e-4, "OK"),
        # Gene 3: Strongly negative wald (wald=-4.5)
        GeneDifferentialExpressionResult("ENSG03", "GENE3", 15.0, -2.25, 0.5, -4.5, 1e-4, 1e-3, "OK"),
        # Gene 4: Invalid stat (NaN) -> MUST BE EXCLUDED
        GeneDifferentialExpressionResult("ENSG04", "GENE4", 5.0, 0.0, None, float("nan"), None, None, "OK"),
        # Gene 5: OUTLIER status -> MUST BE EXCLUDED
        GeneDifferentialExpressionResult("ENSG05", "GENE5", 100.0, 1.0, 0.2, 5.0, 1e-5, 1e-4, "OUTLIER"),
        # Gene 6: Duplicate Entrez 101 with higher absolute wald (wald=3.0 vs Gene 1 wald=1.0)
        GeneDifferentialExpressionResult("ENSG06", "GENE6_DUP", 12.0, 1.5, 0.5, 3.0, 0.01, 0.05, "OK"),
    )

    contrast = DifferentialExpressionContrastResult(
        contrast_id="test_contrast_001",
        study_id="STUDY_01",
        experiment_id="EXP_01",
        cohort_id="COHORT_01",
        relationship_id="REL_01",
        treatment_condition_id="TRT",
        matched_control_condition_ids=("CTL",),
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=("R1", "R2"),
        control_biological_replicate_ids=("R3", "R4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name="compound_a",
        agent_identifier="CID:1234",
        concentration_or_dose=10.0,
        concentration_or_dose_unit="uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 5",
        design_formula="~ condition",
        residual_degrees_of_freedom=2,
        size_factors_used={"S1": 1.0, "S2": 1.0, "S3": 1.0, "S4": 1.0},
        total_input_genes_count=6,
        gene_results=gene_results,
    )
    return contrast, ref


def test_ranking_excludes_invalid_and_preserves_non_deg(
    test_contrast_and_ref: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = test_contrast_and_ref
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    ranked_set, mapping_res, excluded = builder.build_ranking(contrast, deduplicate_entrez=False)

    # Gene 4 (NaN) and Gene 5 (OUTLIER) must be excluded
    excluded_ids = {e.canonical_gene_id for e in excluded}
    assert "ENSG04" in excluded_ids
    assert "ENSG05" in excluded_ids

    # Gene 1 (non-DEG, p=0.4, padj=0.6) MUST BE INCLUDED
    ranked_ids = [g.canonical_gene_id for g in ranked_set.ranked_genes]
    assert "ENSG01" in ranked_ids

    # Ranking values must be ordered strictly decreasingly
    values = [g.ranking_value for g in ranked_set.ranked_genes]
    assert values == sorted(values, reverse=True)


def test_ranking_deduplicates_entrez_by_max_abs_wald(
    test_contrast_and_ref: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = test_contrast_and_ref
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    ranked_set, mapping_res, excluded = builder.build_ranking(contrast, deduplicate_entrez=True)

    # Entrez 101 was shared by ENSG01 (wald=1.0) and ENSG06 (wald=3.0)
    # ENSG06 should be kept (higher |wald|), ENSG01 should be dropped
    ranked_entrez = [g.entrez_gene_id for g in ranked_set.ranked_genes]
    assert ranked_entrez.count("101") == 1

    chosen_gene = next(g for g in ranked_set.ranked_genes if g.entrez_gene_id == "101")
    assert chosen_gene.canonical_gene_id == "ENSG06"
    assert chosen_gene.ranking_value == 3.0

    # Ensure dropped duplicate is recorded in exclusions
    dropped_dup = next(e for e in excluded if e.canonical_gene_id == "ENSG01")
    assert "Duplicate Entrez 101 resolved" in dropped_dup.exclusion_reason


def test_deterministic_ties_handling() -> None:
    records = (
        ReferenceRecord("ENSG_A", "GENE_A", "1", "HGNC:1", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG_B", "GENE_B", "2", "HGNC:2", "protein_coding", (), "1", ""),
    )
    ref = GeneReference("r1", "v1", "GRCh38", "Ensembl", records, "chk1")
    gene_results = (
        GeneDifferentialExpressionResult("ENSG_B", "GENE_B", 10.0, 1.0, 0.5, 2.0, 0.05, 0.1, "OK"),
        GeneDifferentialExpressionResult("ENSG_A", "GENE_A", 10.0, 1.0, 0.5, 2.0, 0.05, 0.1, "OK"),
    )
    contrast = DifferentialExpressionContrastResult(
        contrast_id="c1",
        study_id="s1",
        experiment_id="e1",
        cohort_id="c1",
        relationship_id="r1",
        treatment_condition_id="t",
        matched_control_condition_ids=("c",),
        treatment_sample_ids=("s1",),
        control_sample_ids=("s2",),
        treatment_biological_replicate_ids=("r1",),
        control_biological_replicate_ids=("r2",),
        treatment_replicate_count=1,
        control_replicate_count=1,
        agent_name="compound",
        agent_identifier=None,
        concentration_or_dose=1.0,
        concentration_or_dose_unit="uM",
        exposure_duration=1.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="D1",
        design_formula="~ condition",
        residual_degrees_of_freedom=1,
        size_factors_used={},
        total_input_genes_count=2,
        gene_results=gene_results,
    )
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    ranked_set, _, _ = builder.build_ranking(contrast, deduplicate_entrez=True)
    # Both have ranking_value = 2.0, tie-breaker is canonical_gene_id ascending: ENSG_A before ENSG_B
    assert ranked_set.ranked_genes[0].canonical_gene_id == "ENSG_A"
    assert ranked_set.ranked_genes[1].canonical_gene_id == "ENSG_B"
