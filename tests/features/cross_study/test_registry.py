"""Unit tests for Stage 8A.1 ConditionGeneResponseRegistry.

Verifies:
16. Different VPA doses remain different conditions
17. Two doses from GSE232218 do not become two independent studies
18. Two timepoints from one experiment do not become two independent studies
19. Different compounds cannot be pooled
20. Single-study data return SINGLE_STUDY_ONLY
21. Single-study data do not produce a pooled estimate
22. No heterogeneity statistic (tau^2, I^2) is generated
23. No random-effects meta-analysis runs in Stage 8A.1
26. Row/input ordering does not change outputs
"""

import copy
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)


def test_single_contrast_registry(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200],
        registry_id="reg_vpa_200",
    )

    # 20. Single-study data return SINGLE_STUDY_ONLY
    assert registry.total_records == 8
    assert len(registry.unique_independent_studies) == 1
    assert "GSE232218" in registry.unique_independent_studies

    for r in registry.records:
        assert r.cross_study_status == CrossStudyStatus.SINGLE_STUDY_ONLY

        # 21. No pooled estimate generated
        rec_dict = r.to_dict()
        assert "pooled_log2fc" not in rec_dict
        assert "pooled_se" not in rec_dict

        # 22. No heterogeneity statistic generated
        assert "tau_squared" not in rec_dict
        assert "i_squared" not in rec_dict

    # 23. No meta-analysis executed
    assert registry.provenance.get("meta_analysis_executed") is False


def test_two_doses_remain_distinct_conditions_and_single_study(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
    synthetic_contrast_vpa_400: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200, synthetic_contrast_vpa_400],
        registry_id="reg_vpa_two_doses",
    )

    # 16. Different VPA doses remain different conditions
    assert len(registry.unique_contrasts) == 2
    assert "contrast_GSE232218_VPA_200uM" in registry.unique_contrasts
    assert "contrast_GSE232218_VPA_400uM" in registry.unique_contrasts

    recs_200 = registry.get_records_for_contrast("contrast_GSE232218_VPA_200uM")
    recs_400 = registry.get_records_for_contrast("contrast_GSE232218_VPA_400uM")
    assert len(recs_200) == 8
    assert len(recs_400) == 8

    assert recs_200[0].dose_value == 200.0
    assert recs_400[0].dose_value == 400.0
    assert recs_200[0].treatment_condition_id == "cond_vpa_200"
    assert recs_400[0].treatment_condition_id == "cond_vpa_400"

    # 17. Two doses from GSE232218 do NOT become two independent studies
    assert len(registry.unique_independent_studies) == 1
    assert "GSE232218" in registry.unique_independent_studies

    # Because both contrasts originate from the single study GSE232218, study_count == 1
    assert registry.provenance["independent_study_counts_by_compound"]["CID:3121"] == 1
    for r in registry.records:
        assert r.cross_study_status == CrossStudyStatus.SINGLE_STUDY_ONLY


def test_two_timepoints_from_same_experiment_do_not_become_two_independent_studies(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    # Create second timepoint (24h vs 48h) in the same study
    contrast_24h = copy.deepcopy(synthetic_contrast_vpa_200)
    object.__setattr__(contrast_24h, "contrast_id", "contrast_GSE232218_VPA_24h")
    object.__setattr__(contrast_24h, "exposure_duration", 24.0)
    object.__setattr__(contrast_24h, "treatment_condition_id", "cond_vpa_24h")

    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200, contrast_24h],
        registry_id="reg_timepoints",
    )

    # 18. Two timepoints from one experiment do not become two independent studies
    assert len(registry.unique_contrasts) == 2
    assert len(registry.unique_independent_studies) == 1
    assert registry.provenance["independent_study_counts_by_compound"]["CID:3121"] == 1
    for r in registry.records:
        assert r.cross_study_status == CrossStudyStatus.SINGLE_STUDY_ONLY


def test_different_compounds_cannot_be_pooled(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
    sample_gene_results: tuple[GeneDifferentialExpressionResult, ...],
) -> None:
    # Create contrast for a second compound (Ethanol)
    contrast_etoh = DifferentialExpressionContrastResult(
        contrast_id="contrast_GSE999999_ETOH_50mM",
        study_id="GSE999999",
        experiment_id="exp_02",
        cohort_id="cohort_etoh",
        relationship_id="rel_etoh",
        treatment_condition_id="cond_etoh_50mM",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("S_etoh_1", "S_etoh_2"),
        control_sample_ids=("S_ctrl_1", "S_ctrl_2"),
        treatment_biological_replicate_ids=("R_etoh_1", "R_etoh_2"),
        control_biological_replicate_ids=("R_ctrl_1", "R_ctrl_2"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name="Ethanol",
        agent_identifier="CID:702",
        concentration_or_dose=50.0,
        concentration_or_dose_unit="mM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 15",
        organoid_context_ids=("organoid_ctx_02",),
        biological_source_ids=("WA09",),
        gene_results=sample_gene_results,
    )

    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200, contrast_etoh],
        registry_id="reg_multi_compound",
    )

    # 19. Different compounds cannot be pooled: each compound retains its distinct records
    assert len(registry.unique_compounds) == 2
    assert "CID:3121" in registry.unique_compounds
    assert "CID:702" in registry.unique_compounds

    vpa_recs = registry.get_records_for_compound("CID:3121")
    etoh_recs = registry.get_records_for_compound("CID:702")
    assert len(vpa_recs) == 8
    assert len(etoh_recs) == 8

    # Each compound appears in only 1 independent study -> SINGLE_STUDY_ONLY
    assert registry.provenance["independent_study_counts_by_compound"]["CID:3121"] == 1
    assert registry.provenance["independent_study_counts_by_compound"]["CID:702"] == 1


def test_input_ordering_invariance(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
    synthetic_contrast_vpa_400: DifferentialExpressionContrastResult,
) -> None:
    # 26. Row/input ordering does not change outputs: order contrasts differently
    reg1 = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200, synthetic_contrast_vpa_400],
        registry_id="reg_order1",
    )
    reg2 = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_400, synthetic_contrast_vpa_200],
        registry_id="reg_order2",
    )

    # Both registries should have identical records sorted in canonical order
    assert len(reg1.records) == len(reg2.records)
    for r1, r2 in zip(reg1.records, reg2.records, strict=True):
        assert r1.canonical_ensembl_gene_id == r2.canonical_ensembl_gene_id
        assert r1.contrast_id == r2.contrast_id
        assert r1.log2_fold_change == r2.log2_fold_change
        assert r1.dose_value == r2.dose_value
