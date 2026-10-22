"""Builder functions for converting Stage 7 Differential Expression results into Stage 8A.1 records.

Guarantees 100% numerical and identity preservation:
- Canonical gene ID and symbol
- Exact log2FC, SE, Wald, p-value, and FDR
- Complete chemical, exposure, and biological context
- Explicit missing-value preservation (no replacement with zero)
- Descriptive-only response direction (UP/DOWN/NO_DIRECTION/UNDETERMINED)
"""

import math
from typing import Any

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)


def derive_response_direction(log2_fc: float | None) -> ResponseDirection:
    """Derive purely descriptive direction from the sign of log2FC.

    Non-evaluative: does NOT imply statistical significance or biological toxicity.
    """
    if log2_fc is None or math.isnan(log2_fc):
        return ResponseDirection.UNDETERMINED
    if log2_fc > 0.0:
        return ResponseDirection.UP
    if log2_fc < 0.0:
        return ResponseDirection.DOWN
    return ResponseDirection.NO_DIRECTION


def derive_estimate_status(gene_result: GeneDifferentialExpressionResult) -> tuple[EstimateStatus, str | None]:
    """Classify the estimation status and reason for missing/invalid statistics."""
    # 1. Zero counts check
    if gene_result.result_status == "ZERO_COUNTS" or gene_result.base_mean == 0.0:
        return EstimateStatus.ZERO_COUNTS, "baseMean == 0 across all samples in contrast"

    # 2. Outlier check (DESeq2 sets pvalue to NA for Cook's distance outliers)
    if gene_result.result_status == "OUTLIER" or (
        gene_result.p_value is None
        and gene_result.log2_fold_change is not None
        and gene_result.base_mean > 0.0
    ):
        return (
            EstimateStatus.COOKS_OUTLIER,
            "Cook's distance outlier; p-value set to NA by DESeq2",
        )

    # 3. Fit error check
    if gene_result.result_status == "FIT_ERROR":
        return EstimateStatus.FIT_ERROR, gene_result.result_notes or "Model fitting failure in DESeq2"

    # 4. Valid estimate check
    if (
        gene_result.log2_fold_change is not None
        and gene_result.lfc_standard_error is not None
        and gene_result.wald_statistic is not None
        and gene_result.p_value is not None
    ):
        if not (math.isnan(gene_result.wald_statistic) or math.isinf(gene_result.wald_statistic)):
            return EstimateStatus.OK, None
        return EstimateStatus.INVALID_STATISTIC, "Wald statistic is NaN or infinite"

    return EstimateStatus.INVALID_STATISTIC, gene_result.result_notes or "Incomplete or missing statistical estimates"


def _parse_dose_value(dose_input: float | str | None) -> float | None:
    """Safely extract numeric dose value without fabricating numbers."""
    if dose_input is None:
        return None
    if isinstance(dose_input, (int, float)):
        return float(dose_input)
    try:
        # If it's a string like "200.0" or "200"
        return float(dose_input)
    except (ValueError, TypeError):
        return None


def _parse_duration_value(dur_input: float | str | None) -> float | None:
    """Safely extract numeric exposure duration without fabricating numbers."""
    if dur_input is None:
        return None
    if isinstance(dur_input, (int, float)):
        return float(dur_input)
    try:
        return float(dur_input)
    except (ValueError, TypeError):
        return None


def normalize_independent_study_id(study_id: str, study_family_id: str | None = None) -> str:
    """Derive standardized independent study identifier.

    Multiple conditions/contrasts from the same study accession/family share this identifier.
    """
    if study_family_id and study_family_id.strip():
        return study_family_id.strip()
    clean = study_id.strip()
    if clean.startswith("study_"):
        clean = clean[len("study_") :]
    return clean


def build_condition_gene_response_record(
    gene_result: GeneDifferentialExpressionResult,
    contrast: DifferentialExpressionContrastResult,
    study_family_id: str | None = None,
    model_type: str | None = None,
    cross_study_status: CrossStudyStatus = CrossStudyStatus.SINGLE_STUDY_ONLY,
    registry_version: str = "1.0.0",
) -> ConditionGeneResponseRecord:
    """Construct an immutable ConditionGeneResponseRecord from a Stage 7 gene result and contrast."""
    ind_study_id = normalize_independent_study_id(contrast.study_id, study_family_id)
    fam_id = study_family_id.strip() if study_family_id and study_family_id.strip() else contrast.study_id

    # Compound context
    comp_id = contrast.agent_identifier or contrast.agent_name or "UNKNOWN_COMPOUND"
    comp_name = contrast.agent_name or contrast.agent_identifier or "UNKNOWN_COMPOUND"

    # Exposure context
    dose_val = _parse_dose_value(contrast.concentration_or_dose)
    dose_unit = contrast.concentration_or_dose_unit
    dose_orig = str(contrast.concentration_or_dose) if contrast.concentration_or_dose is not None else None

    dur_val = _parse_duration_value(contrast.exposure_duration)
    dur_unit = contrast.exposure_duration_unit

    dev_stage = contrast.developmental_age_or_stage_at_exposure or contrast.treatment_collection_age_or_stage

    # Direction and estimate status
    direction = derive_response_direction(gene_result.log2_fold_change)
    est_status, missing_reason = derive_estimate_status(gene_result)

    return ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=gene_result.canonical_gene_id,
        gene_symbol=gene_result.approved_symbol,
        study_id=contrast.study_id,
        study_family_id=fam_id,
        independent_study_id=ind_study_id,
        experiment_id=contrast.experiment_id,
        contrast_id=contrast.contrast_id,
        treatment_condition_id=contrast.treatment_condition_id,
        matched_control_condition_ids=contrast.matched_control_condition_ids,
        compound_id=comp_id,
        compound_name=comp_name,
        dose_value=dose_val,
        dose_unit=dose_unit,
        dose_original=dose_orig,
        exposure_duration=dur_val,
        exposure_duration_unit=dur_unit,
        developmental_stage=dev_stage,
        model_type=model_type or (contrast.organoid_context_ids[0] if contrast.organoid_context_ids else None),
        biological_source_ids=contrast.biological_source_ids,
        organoid_context_ids=contrast.organoid_context_ids,
        treatment_replicate_count=contrast.treatment_replicate_count,
        control_replicate_count=contrast.control_replicate_count,
        treatment_sample_ids=contrast.treatment_sample_ids,
        control_sample_ids=contrast.control_sample_ids,
        treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
        control_biological_replicate_ids=contrast.control_biological_replicate_ids,
        base_mean=gene_result.base_mean,
        log2_fold_change=gene_result.log2_fold_change,
        log2_fold_change_se=gene_result.lfc_standard_error,
        wald_statistic=gene_result.wald_statistic,
        p_value=gene_result.p_value,
        adjusted_p_value=gene_result.adjusted_p_value_bh,
        response_direction=direction,
        estimate_status=est_status,
        missing_reason=missing_reason,
        cross_study_status=cross_study_status,
        source_stage="Stage 7B",
        source_contrast_id=contrast.contrast_id,
        de_version=contrast.de_version,
        registry_version=registry_version,
        input_hashes=dict(contrast.input_hashes),
    )


def build_records_for_contrast(
    contrast: DifferentialExpressionContrastResult,
    study_family_id: str | None = None,
    model_type: str | None = None,
    cross_study_status: CrossStudyStatus = CrossStudyStatus.SINGLE_STUDY_ONLY,
    registry_version: str = "1.0.0",
) -> tuple[ConditionGeneResponseRecord, ...]:
    """Convert all gene results in a Stage 7 contrast into ConditionGeneResponseRecord objects."""
    records = []
    for g in contrast.gene_results:
        rec = build_condition_gene_response_record(
            gene_result=g,
            contrast=contrast,
            study_family_id=study_family_id,
            model_type=model_type,
            cross_study_status=cross_study_status,
            registry_version=registry_version,
        )
        records.append(rec)
    return tuple(records)
