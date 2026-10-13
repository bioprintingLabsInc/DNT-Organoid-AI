"""Tests verifying strict prevention of DNT outcome label leakage and circularity in Stage 8C."""

from pathlib import Path
import re
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.pathway.mapping import EnsemblToEntrezAdapter
from src.features.pathway.ranking import GeneRankingBuilder
from src.gene_harmonization.reference import GeneReference, ReferenceRecord


def test_no_dnt_labels_in_stage8c_source_code() -> None:
    """Audit all source files in src/features/pathway/ for forbidden outcome-label tokens."""
    pathway_src_dir = Path("src/features/pathway")
    assert pathway_src_dir.is_dir()

    forbidden_patterns = [
        re.compile(r"\bdnt_positive\b", re.IGNORECASE),
        re.compile(r"\bdnt_negative\b", re.IGNORECASE),
        re.compile(r"\btoxicant_class\b", re.IGNORECASE),
        re.compile(r"\bclassifier\b", re.IGNORECASE),
        re.compile(r"\belastic_net\b", re.IGNORECASE),
        re.compile(r"\brandom_forest\b", re.IGNORECASE),
        re.compile(r"\bbglmm\b", re.IGNORECASE),
        re.compile(r"\bprediction_label\b", re.IGNORECASE),
    ]

    for py_file in pathway_src_dir.glob("*.py"):
        content = py_file.read_text(encoding="utf-8")
        for pat in forbidden_patterns:
            matches = pat.findall(content)
            assert not matches, f"Found forbidden leakage pattern '{pat.pattern}' in {py_file}"


def test_ranking_contract_does_not_accept_or_require_outcome_labels() -> None:
    """Ensure GeneRankingBuilder operates purely on statistical and experimental metadata."""
    records = (
        ReferenceRecord("ENSG01", "G1", "101", "HGNC:1", "protein_coding", (), "1", ""),
    )
    ref = GeneReference("r1", "v1", "GRCh38", "Ensembl", records, "chk1")
    gene_results = (
        GeneDifferentialExpressionResult("ENSG01", "G1", 10.0, 1.0, 0.2, 5.0, 1e-5, 1e-4, "OK"),
    )
    contrast = DifferentialExpressionContrastResult(
        contrast_id="c_unlabeled",
        study_id="STUDY_X",
        experiment_id="EXP_X",
        cohort_id="COHORT_X",
        relationship_id="REL_X",
        treatment_condition_id="TRT_X",
        matched_control_condition_ids=("CTL_X",),
        treatment_sample_ids=("S1",),
        control_sample_ids=("S2",),
        treatment_biological_replicate_ids=("R1",),
        control_biological_replicate_ids=("R2",),
        treatment_replicate_count=1,
        control_replicate_count=1,
        agent_name="mystery_compound_without_class_label",
        agent_identifier=None,
        concentration_or_dose=5.0,
        concentration_or_dose_unit="uM",
        exposure_duration=12.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 3",
        design_formula="~ condition",
        residual_degrees_of_freedom=1,
        size_factors_used={},
        total_input_genes_count=1,
        gene_results=gene_results,
    )

    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)
    ranked_set, mapping_res, excluded = builder.build_ranking(contrast)

    # Ranking successfully constructed without any knowledge of compound outcome
    assert len(ranked_set.ranked_genes) == 1
    assert ranked_set.ranked_genes[0].canonical_gene_id == "ENSG01"
    assert ranked_set.ranked_genes[0].ranking_value == 5.0
