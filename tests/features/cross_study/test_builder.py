"""Unit tests for Stage 8A.1 builder functions.

Verifies:
1. Canonical gene identity preservation
2. log2FC, SE, Wald, p-value, and FDR numerical preservation
3. Direction logic (positive -> UP, negative -> DOWN, zero -> NO_DIRECTION, None/NaN -> UNDETERMINED)
4. Absence of magnitude cutoffs (no minimum |log2FC| threshold)
5. Absence of significance cutoffs (genes with FDR > 0.05 remain preserved)
6. Non-Entrez gene preservation
7. Missing-data traceability and non-fabrication of missing metadata
8. Correct treatment vs matched control assignment
"""

import math
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.builder import (
    build_condition_gene_response_record,
    build_records_for_contrast,
    derive_estimate_status,
    derive_response_direction,
)
from src.features.cross_study.models import (
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)


def test_derive_response_direction() -> None:
    # 7. positive log2FC -> UP
    assert derive_response_direction(1.5) == ResponseDirection.UP
    assert derive_response_direction(0.001) == ResponseDirection.UP

    # 8. negative log2FC -> DOWN
    assert derive_response_direction(-1.5) == ResponseDirection.DOWN
    assert derive_response_direction(-0.0001) == ResponseDirection.DOWN

    # 9. zero log2FC -> NO_DIRECTION
    assert derive_response_direction(0.0) == ResponseDirection.NO_DIRECTION

    # None and NaN -> UNDETERMINED
    assert derive_response_direction(None) == ResponseDirection.UNDETERMINED
    assert derive_response_direction(float("nan")) == ResponseDirection.UNDETERMINED


def test_derive_estimate_status() -> None:
    # Valid estimate
    g_ok = GeneDifferentialExpressionResult("G1", "SYM1", 10.0, 1.0, 0.2, 5.0, 1e-4, 1e-3, "OK")
    st, reason = derive_estimate_status(g_ok)
    assert st == EstimateStatus.OK
    assert reason is None

    # Zero counts
    g_zero = GeneDifferentialExpressionResult("G2", "SYM2", 0.0, None, None, None, None, None, "ZERO_COUNTS")
    st_z, reason_z = derive_estimate_status(g_zero)
    assert st_z == EstimateStatus.ZERO_COUNTS
    assert "baseMean == 0" in (reason_z or "")

    # Cook's outlier (p-value is None, log2FC exists, baseMean > 0)
    g_out = GeneDifferentialExpressionResult("G3", "SYM3", 100.0, 2.0, 0.5, 4.0, None, None, "OUTLIER")
    st_o, reason_o = derive_estimate_status(g_out)
    assert st_o == EstimateStatus.COOKS_OUTLIER
    assert "Cook's distance" in (reason_o or "")

    # Fit error
    g_err = GeneDifferentialExpressionResult("G4", "SYM4", 50.0, None, None, None, None, None, "FIT_ERROR")
    st_e, reason_e = derive_estimate_status(g_err)
    assert st_e == EstimateStatus.FIT_ERROR

    # NaN wald statistic
    g_nan = GeneDifferentialExpressionResult("G5", "SYM5", 50.0, 1.0, 0.5, float("nan"), 0.05, 0.1, "OK")
    st_n, reason_n = derive_estimate_status(g_nan)
    assert st_n == EstimateStatus.INVALID_STATISTIC


def test_preservation_of_all_genes_and_numbers(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    records = build_records_for_contrast(synthetic_contrast_vpa_200)

    # 1. Canonical gene identity preserved for all 8 input genes
    assert len(records) == len(synthetic_contrast_vpa_200.gene_results)
    record_map = {r.canonical_ensembl_gene_id: r for r in records}
    contrast_map = {g.canonical_gene_id: g for g in synthetic_contrast_vpa_200.gene_results}

    for gid, g in contrast_map.items():
        assert gid in record_map
        r = record_map[gid]

        # 1. Canonical ID & Symbol
        assert r.canonical_ensembl_gene_id == g.canonical_gene_id
        assert r.gene_symbol == g.approved_symbol

        # 2. log2FC is unchanged
        assert r.log2_fold_change == g.log2_fold_change

        # 3. SE is unchanged
        assert r.log2_fold_change_se == g.lfc_standard_error

        # 4. Wald statistic is unchanged
        assert r.wald_statistic == g.wald_statistic

        # 5. p-value is unchanged
        assert r.p_value == g.p_value

        # 6. FDR / adjusted p-value is unchanged
        assert r.adjusted_p_value == g.adjusted_p_value_bh

        # baseMean is unchanged
        assert r.base_mean == g.base_mean


def test_no_arbitrary_magnitude_or_significance_cutoffs(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    records = build_records_for_contrast(synthetic_contrast_vpa_200)
    record_map = {r.canonical_ensembl_gene_id: r for r in records}

    # 10. No absolute log2FC threshold: small fold change gene remains present
    r_small = record_map["ENSG00000000004"]
    assert abs(r_small.log2_fold_change or 0.0) < 0.2
    assert r_small.estimate_status == EstimateStatus.OK

    # 11. Non-significant genes remain present: FDR > 0.05
    r_nonsig = record_map["ENSG00000000002"]
    assert r_nonsig.adjusted_p_value == 0.652
    assert r_nonsig.estimate_status == EstimateStatus.OK

    # 12. Non-Entrez gene is allowed: lncRNA ENSG00000284742
    r_lnc = record_map["ENSG00000284742"]
    assert r_lnc.canonical_ensembl_gene_id == "ENSG00000284742"
    assert r_lnc.log2_fold_change == 1.20
    assert r_lnc.estimate_status == EstimateStatus.OK


def test_missing_statistics_and_metadata_not_fabricated(
    sample_gene_results: tuple[GeneDifferentialExpressionResult, ...],
) -> None:
    # Construct contrast with some missing metadata fields
    contrast_sparse = DifferentialExpressionContrastResult(
        contrast_id="contrast_sparse_01",
        study_id="STUDY_SPARSE",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="cond_trt",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=(),
        control_biological_replicate_ids=(),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name=None,
        agent_identifier="CID:9999",
        concentration_or_dose=None,  # Missing dose
        concentration_or_dose_unit=None,
        exposure_duration=None,  # Missing duration
        exposure_duration_unit=None,
        developmental_age_or_stage_at_exposure=None,  # Missing stage
        organoid_context_ids=(),
        biological_source_ids=(),
        gene_results=sample_gene_results,
    )

    records = build_records_for_contrast(contrast_sparse)
    r_first = records[0]

    # 14. Missing metadata is not fabricated (remains None, not replaced with 0 or empty string)
    assert r_first.dose_value is None
    assert r_first.dose_unit is None
    assert r_first.dose_original is None
    assert r_first.exposure_duration is None
    assert r_first.exposure_duration_unit is None
    assert r_first.developmental_stage is None

    # 13. Missing statistics remain traceable (e.g. zero count gene)
    record_map = {r.canonical_ensembl_gene_id: r for r in records}
    r_zero = record_map["ENSG00000000008"]
    assert r_zero.estimate_status == EstimateStatus.ZERO_COUNTS
    assert r_zero.log2_fold_change is None
    assert r_zero.missing_reason is not None

    # 15. Treatment and control condition identity preserved
    assert r_first.treatment_condition_id == "cond_trt"
    assert r_first.matched_control_condition_ids == ("cond_ctrl",)
    assert r_first.contrast_id == "contrast_sparse_01"
