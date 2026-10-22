"""Pytest fixtures for Stage 8A.1 Cross-Study Gene Response tests."""

import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)


@pytest.fixture
def sample_gene_results() -> tuple[GeneDifferentialExpressionResult, ...]:
    """Synthetic gene differential expression results covering all estimation states and directions."""
    return (
        # 1. UP gene, statistically significant (FDR < 0.05)
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000001",
            approved_symbol="GENE_UP_SIG",
            base_mean=120.5,
            log2_fold_change=1.85,
            lfc_standard_error=0.25,
            wald_statistic=7.40,
            p_value=1.36e-13,
            adjusted_p_value_bh=2.5e-11,
            result_status="OK",
        ),
        # 2. UP gene, NOT significant (FDR > 0.05) -> MUST REMAIN IN REGISTRY
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000002",
            approved_symbol="GENE_UP_NONSIG",
            base_mean=45.2,
            log2_fold_change=0.35,
            lfc_standard_error=0.30,
            wald_statistic=1.167,
            p_value=0.243,
            adjusted_p_value_bh=0.652,
            result_status="OK",
        ),
        # 3. DOWN gene, statistically significant
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000003",
            approved_symbol="GENE_DOWN_SIG",
            base_mean=88.0,
            log2_fold_change=-2.10,
            lfc_standard_error=0.35,
            wald_statistic=-6.00,
            p_value=1.97e-9,
            adjusted_p_value_bh=1.5e-7,
            result_status="OK",
        ),
        # 4. DOWN gene, NOT significant
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000004",
            approved_symbol="GENE_DOWN_NONSIG",
            base_mean=32.1,
            log2_fold_change=-0.15,
            lfc_standard_error=0.20,
            wald_statistic=-0.75,
            p_value=0.453,
            adjusted_p_value_bh=0.812,
            result_status="OK",
        ),
        # 5. NO_DIRECTION gene (log2FC == 0.0)
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000005",
            approved_symbol="GENE_ZERO_FC",
            base_mean=50.0,
            log2_fold_change=0.0,
            lfc_standard_error=0.10,
            wald_statistic=0.0,
            p_value=1.0,
            adjusted_p_value_bh=1.0,
            result_status="OK",
        ),
        # 6. Non-Entrez valid gene (lncRNA without Entrez) -> MUST BE ALLOWED AND PRESERVED
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000284742",
            approved_symbol="LINC0001",
            base_mean=15.0,
            log2_fold_change=1.20,
            lfc_standard_error=0.40,
            wald_statistic=3.00,
            p_value=0.0027,
            adjusted_p_value_bh=0.015,
            result_status="OK",
        ),
        # 7. Cook's distance outlier (p-value and padj are NA)
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000007",
            approved_symbol="GENE_OUTLIER",
            base_mean=500.0,
            log2_fold_change=3.5,
            lfc_standard_error=0.8,
            wald_statistic=4.375,
            p_value=None,
            adjusted_p_value_bh=None,
            result_status="OUTLIER",
            result_notes="Cook's distance outlier",
        ),
        # 8. All-zero count gene (baseMean == 0)
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000008",
            approved_symbol="GENE_ALL_ZERO",
            base_mean=0.0,
            log2_fold_change=None,
            lfc_standard_error=None,
            wald_statistic=None,
            p_value=None,
            adjusted_p_value_bh=None,
            result_status="ZERO_COUNTS",
            result_notes="baseMean == 0",
        ),
    )


@pytest.fixture
def synthetic_contrast_vpa_200(
    sample_gene_results: tuple[GeneDifferentialExpressionResult, ...],
) -> DifferentialExpressionContrastResult:
    """Synthetic contrast for GSE232218 VPA 200 uM vs Vehicle."""
    return DifferentialExpressionContrastResult(
        contrast_id="contrast_GSE232218_VPA_200uM",
        study_id="GSE232218",
        experiment_id="exp_01",
        cohort_id="cohort_vpa_200",
        relationship_id="rel_01",
        treatment_condition_id="cond_vpa_200",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("S_trt_1", "S_trt_2", "S_trt_3", "S_trt_4"),
        control_sample_ids=("S_ctrl_1", "S_ctrl_2", "S_ctrl_3", "S_ctrl_4"),
        treatment_biological_replicate_ids=("R_trt_1", "R_trt_2", "R_trt_3", "R_trt_4"),
        control_biological_replicate_ids=("R_ctrl_1", "R_ctrl_2", "R_ctrl_3", "R_ctrl_4"),
        treatment_replicate_count=4,
        control_replicate_count=4,
        agent_name="Valproic Acid",
        agent_identifier="CID:3121",
        vehicle="water",
        concentration_or_dose=200.0,
        concentration_or_dose_unit="uM",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 6",
        organoid_context_ids=("organoid_ctx_01",),
        biological_source_ids=("WA09",),
        design_formula="~ condition",
        size_factors_used={f"S_trt_{i}": 1.0 for i in range(1, 5)} | {f"S_ctrl_{i}": 1.0 for i in range(1, 5)},
        total_input_genes_count=len(sample_gene_results),
        gene_results=sample_gene_results,
    )


@pytest.fixture
def synthetic_contrast_vpa_400(
    sample_gene_results: tuple[GeneDifferentialExpressionResult, ...],
) -> DifferentialExpressionContrastResult:
    """Second dose condition in the same study: GSE232218 VPA 400 uM vs Vehicle."""
    # Modulate fold changes slightly to simulate higher dose
    modified_results = []
    for g in sample_gene_results:
        lfc = g.log2_fold_change * 1.5 if g.log2_fold_change is not None else None
        wald = g.wald_statistic * 1.4 if g.wald_statistic is not None else None
        modified_results.append(
            GeneDifferentialExpressionResult(
                canonical_gene_id=g.canonical_gene_id,
                approved_symbol=g.approved_symbol,
                base_mean=g.base_mean,
                log2_fold_change=lfc,
                lfc_standard_error=g.lfc_standard_error,
                wald_statistic=wald,
                p_value=g.p_value,
                adjusted_p_value_bh=g.adjusted_p_value_bh,
                result_status=g.result_status,
                result_notes=g.result_notes,
            )
        )

    return DifferentialExpressionContrastResult(
        contrast_id="contrast_GSE232218_VPA_400uM",
        study_id="GSE232218",
        experiment_id="exp_01",
        cohort_id="cohort_vpa_400",
        relationship_id="rel_02",
        treatment_condition_id="cond_vpa_400",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("S_trt_5", "S_trt_6", "S_trt_7", "S_trt_8"),
        control_sample_ids=("S_ctrl_1", "S_ctrl_2", "S_ctrl_3", "S_ctrl_4"),
        treatment_biological_replicate_ids=("R_trt_5", "R_trt_6", "R_trt_7", "R_trt_8"),
        control_biological_replicate_ids=("R_ctrl_1", "R_ctrl_2", "R_ctrl_3", "R_ctrl_4"),
        treatment_replicate_count=4,
        control_replicate_count=4,
        agent_name="Valproic Acid",
        agent_identifier="CID:3121",
        vehicle="water",
        concentration_or_dose=400.0,
        concentration_or_dose_unit="uM",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 6",
        organoid_context_ids=("organoid_ctx_01",),
        biological_source_ids=("WA09",),
        design_formula="~ condition",
        size_factors_used={f"S_trt_{i}": 1.0 for i in range(5, 9)} | {f"S_ctrl_{i}": 1.0 for i in range(1, 5)},
        total_input_genes_count=len(modified_results),
        gene_results=tuple(modified_results),
    )
