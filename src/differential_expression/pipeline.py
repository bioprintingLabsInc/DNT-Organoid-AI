"""Pipeline orchestration for Step 7B Bulk Differential Expression v1."""

from collections.abc import Sequence
import math
from typing import Any

from src.gene_harmonization.models import GeneMappingStatus, HarmonizedDataset
from src.normalization.models import NormalizedCountMatrix, NormalizedDataset
from src.qc.matrix import CountMatrix
from src.response_builder.models import (
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
    ordered_findings,
    status_for,
)

from .design import DesignResolver
from .models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    ExcludedGeneAudit,
    GeneDifferentialExpressionResult,
)
from .r_bridge import execute_deseq2_contrast_r


def run_contrast_differential_expression(
    contrast: MolecularResponseContrast,
    normalized_dataset: NormalizedDataset,
    harmonized_dataset: HarmonizedDataset,
    sample_metadata: dict[str, Any] | Sequence[dict[str, Any]] | None = None,
    r_binary_path: str = "Rscript",
    strict_version_check: bool = True,
    r_timeout_seconds: int = 180,
    config_version: str = "1.0.0",
    de_version: str = "1.0.0",
) -> DifferentialExpressionContrastResult:
    """Execute differential expression for a single Step 7A MolecularResponseContrast."""
    findings: list[Finding] = []

    # 1. Gating on Step 7A contrast qualification
    if contrast.status != ContrastStatus.ELIGIBLE or not contrast.is_inferentially_eligible:
        findings.extend(contrast.findings)
        if not contrast.is_inferentially_eligible and contrast.status == ContrastStatus.ELIGIBLE:
            res_status = ContrastStatus.BLOCKED
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="insufficient_biological_replicates",
                    entity_type="contrast",
                    entity_id=contrast.contrast_id,
                    path="biological_replicate_count",
                    message=(
                        f"Contrast '{contrast.contrast_id}' has insufficient biological replicates "
                        f"(treatment={contrast.treatment_replicate_count}, control={contrast.control_replicate_count}); "
                        "minimum 2 required in each arm for differential expression inference."
                    ),
                )
            )
        else:
            res_status = contrast.status
            findings.append(
                Finding(
                    severity=Severity.INFO,
                    rule_id="contrast_ineligible_for_de",
                    entity_type="contrast",
                    entity_id=contrast.contrast_id,
                    path="status",
                    message=(
                        f"Contrast '{contrast.contrast_id}' has status {contrast.status} "
                        f"(is_inferentially_eligible={contrast.is_inferentially_eligible}); "
                        "DESeq2 execution skipped."
                    ),
                )
            )
        return DifferentialExpressionContrastResult(
            contrast_id=contrast.contrast_id,
            study_id=contrast.study_id,
            experiment_id=contrast.experiment_id,
            cohort_id=contrast.cohort_id,
            relationship_id=contrast.relationship_id,
            treatment_condition_id=contrast.treatment_condition_id,
            matched_control_condition_ids=contrast.matched_control_condition_ids,
            treatment_sample_ids=contrast.treatment_sample_ids,
            control_sample_ids=contrast.control_sample_ids,
            treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
            control_biological_replicate_ids=contrast.control_biological_replicate_ids,
            treatment_replicate_count=contrast.treatment_replicate_count,
            control_replicate_count=contrast.control_replicate_count,
            treatment_exposures=contrast.treatment_exposures,
            sample_level_exposure_deviations=contrast.sample_level_exposure_deviations,
            agent_name=contrast.agent_name,
            agent_identifier=contrast.agent_identifier,
            vehicle=contrast.vehicle,
            concentration_or_dose=contrast.concentration_or_dose,
            concentration_or_dose_unit=contrast.concentration_or_dose_unit,
            exposure_start_time_or_stage=contrast.exposure_start_time_or_stage,
            developmental_age_or_stage_at_exposure=contrast.developmental_age_or_stage_at_exposure,
            exposure_duration=contrast.exposure_duration,
            exposure_duration_unit=contrast.exposure_duration_unit,
            washout_or_recovery_duration=contrast.washout_or_recovery_duration,
            washout_or_recovery_duration_unit=contrast.washout_or_recovery_duration_unit,
            treatment_collection_age_or_stage=contrast.treatment_collection_age_or_stage,
            control_collection_age_or_stage=contrast.control_collection_age_or_stage,
            treatment_collection_age_normalized=contrast.treatment_collection_age_normalized,
            control_collection_age_normalized=contrast.control_collection_age_normalized,
            biological_source_ids=contrast.biological_source_ids,
            organoid_context_ids=contrast.organoid_context_ids,
            design_formula="",
            design_rationale=f"Skipped upstream: contrast status is {res_status} (is_inferentially_eligible={contrast.is_inferentially_eligible}).",
            status=res_status,
            findings=ordered_findings(findings),
            configuration_version=config_version,
            de_version=de_version,
        )

    # 2. Find corresponding cohort in NormalizedDataset
    cohort_res = None
    for cr in normalized_dataset.cohort_results:
        if cr.cohort.cohort_id == contrast.cohort_id:
            cohort_res = cr
            break

    if cohort_res is None:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="cohort_not_found",
                entity_type="cohort",
                entity_id=contrast.cohort_id,
                path="cohort_id",
                message=f"Normalization cohort '{contrast.cohort_id}' not found in normalized dataset.",
            )
        )
        return _make_blocked_result(contrast, findings, config_version, de_version, "Cohort not found in normalized dataset.")

    # 3. Validate Raw Integer Count Matrix
    raw_matrix = cohort_res.raw_count_matrix
    if isinstance(raw_matrix, NormalizedCountMatrix) or not isinstance(raw_matrix, CountMatrix):
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="invalid_count_matrix_type",
                entity_type="count_matrix",
                entity_id=contrast.cohort_id,
                path="raw_count_matrix",
                message="Differential expression requires raw integer CountMatrix, not normalized floating-point counts.",
            )
        )
        return _make_blocked_result(contrast, findings, config_version, de_version, "Invalid count matrix type.")

    # Verify counts are integers >= 0
    try:
        for col in raw_matrix.columns:
            for val in col:
                if not isinstance(val, int) or val < 0:
                    raise ValueError(f"Non-integer or negative count: {val}")
    except (ValueError, TypeError) as e:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="invalid_raw_counts",
                entity_type="count_matrix",
                entity_id=contrast.cohort_id,
                path="raw_count_matrix",
                message=f"Raw count matrix contains non-integer or negative values: {e}",
            )
        )
        return _make_blocked_result(contrast, findings, config_version, de_version, "Raw counts must be non-negative integers.")

    # 4. Size Factor Validation
    all_contrast_samples = tuple(sorted(set(contrast.control_sample_ids + contrast.treatment_sample_ids)))
    sf_by_sample: dict[str, float] = {}
    for sf_item in cohort_res.size_factors:
        sf_by_sample[sf_item.sample_id] = sf_item.size_factor

    missing_sf = [s for s in all_contrast_samples if s not in sf_by_sample]
    if missing_sf:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="missing_size_factors",
                entity_type="size_factor",
                entity_id=contrast.cohort_id,
                path="size_factors",
                message=f"Missing locked size factors for samples: {missing_sf}.",
            )
        )
        return _make_blocked_result(contrast, findings, config_version, de_version, "Missing locked size factors.")

    invalid_sf = [
        s for s in all_contrast_samples
        if sf_by_sample[s] is None or not math.isfinite(sf_by_sample[s]) or sf_by_sample[s] <= 0.0
    ]
    if invalid_sf:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="invalid_size_factors",
                entity_type="size_factor",
                entity_id=contrast.cohort_id,
                path="size_factors",
                message=f"Size factors must be finite and strictly positive; invalid for samples: {invalid_sf}.",
            )
        )
        return _make_blocked_result(contrast, findings, config_version, de_version, "Invalid size factor values.")

    size_factors_used = {s: sf_by_sample[s] for s in all_contrast_samples}

    # 5. Gene Harmonization Mapping & Gene-Universe Conservation Audit
    total_input_genes = len(raw_matrix.gene_ids)
    colliding_canonical_ids = {c.canonical_gene_id for c in harmonized_dataset.collisions}

    # Count occurrences of uniquely mapped canonical gene IDs to catch any collisions
    canonical_counts: dict[str, int] = {}
    for g in harmonized_dataset.genes:
        if g.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED and g.canonical_gene_id:
            canonical_counts[g.canonical_gene_id] = canonical_counts.get(g.canonical_gene_id, 0) + 1

    harm_by_idx = {g.source_index: g for g in harmonized_dataset.genes}
    harm_by_orig = {g.original_gene_id: g for g in harmonized_dataset.genes}

    raw_canonical_counts: dict[str, int] = {}
    for idx, orig_id in enumerate(raw_matrix.gene_ids):
        harm_gene = None
        if idx in harm_by_idx and harm_by_idx[idx].original_gene_id == orig_id:
            harm_gene = harm_by_idx[idx]
        elif orig_id in harm_by_orig:
            harm_gene = harm_by_orig[orig_id]
        if harm_gene and harm_gene.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED and harm_gene.canonical_gene_id:
            raw_canonical_counts[harm_gene.canonical_gene_id] = raw_canonical_counts.get(harm_gene.canonical_gene_id, 0) + 1

    all_colliding_canonicals = (
        colliding_canonical_ids
        | {cid for cid, count in canonical_counts.items() if count > 1}
        | {cid for cid, count in raw_canonical_counts.items() if count > 1}
    )

    valid_entries: list[tuple[str, str, int]] = []
    symbol_map: dict[str, str | None] = {}
    excluded_gene_audits: list[ExcludedGeneAudit] = []

    for idx, orig_id in enumerate(raw_matrix.gene_ids):
        harm_gene = None
        if idx in harm_by_idx and harm_by_idx[idx].original_gene_id == orig_id:
            harm_gene = harm_by_idx[idx]
        elif orig_id in harm_by_orig:
            harm_gene = harm_by_orig[orig_id]

        if harm_gene is None:
            excluded_gene_audits.append(
                ExcludedGeneAudit(
                    original_gene_id=orig_id,
                    source_index=idx,
                    mapping_status="UNMAPPED",
                    canonical_gene_id=None,
                    exclusion_reason="Gene identifier not found in harmonized dataset.",
                )
            )
        elif harm_gene.mapping_status == GeneMappingStatus.UNMAPPED:
            excluded_gene_audits.append(
                ExcludedGeneAudit(
                    original_gene_id=orig_id,
                    source_index=idx,
                    mapping_status="UNMAPPED",
                    canonical_gene_id=None,
                    exclusion_reason=harm_gene.mapping_reason or "Unmapped gene identifier.",
                )
            )
        elif harm_gene.mapping_status == GeneMappingStatus.AMBIGUOUS:
            excluded_gene_audits.append(
                ExcludedGeneAudit(
                    original_gene_id=orig_id,
                    source_index=idx,
                    mapping_status="AMBIGUOUS",
                    canonical_gene_id=harm_gene.canonical_gene_id,
                    exclusion_reason=harm_gene.mapping_reason or "Ambiguous gene identifier mapping.",
                )
            )
        elif harm_gene.mapping_status == GeneMappingStatus.INVALID:
            excluded_gene_audits.append(
                ExcludedGeneAudit(
                    original_gene_id=orig_id,
                    source_index=idx,
                    mapping_status="INVALID",
                    canonical_gene_id=None,
                    exclusion_reason=harm_gene.mapping_reason or "Invalid gene identifier.",
                )
            )
        elif harm_gene.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED:
            can_id = harm_gene.canonical_gene_id
            if not can_id:
                excluded_gene_audits.append(
                    ExcludedGeneAudit(
                        original_gene_id=orig_id,
                        source_index=idx,
                        mapping_status="INVALID",
                        canonical_gene_id=None,
                        exclusion_reason="Uniquely mapped status but canonical gene ID is null.",
                    )
                )
            elif can_id in all_colliding_canonicals:
                excluded_gene_audits.append(
                    ExcludedGeneAudit(
                        original_gene_id=orig_id,
                        source_index=idx,
                        mapping_status="COLLISION",
                        canonical_gene_id=can_id,
                        exclusion_reason=f"Excluded due to canonical collision: multiple source identifiers map to '{can_id}'.",
                    )
                )
            else:
                valid_entries.append((can_id, orig_id, idx))
                symbol_map[can_id] = harm_gene.approved_symbol
        else:
            excluded_gene_audits.append(
                ExcludedGeneAudit(
                    original_gene_id=orig_id,
                    source_index=idx,
                    mapping_status=str(harm_gene.mapping_status),
                    canonical_gene_id=harm_gene.canonical_gene_id,
                    exclusion_reason=harm_gene.mapping_reason or f"Excluded due to unhandled mapping status: {harm_gene.mapping_status}",
                )
            )

    if excluded_gene_audits:
        findings.append(
            Finding(
                severity=Severity.INFO,
                rule_id="gene_universe_harmonization_exclusions",
                entity_type="gene_harmonization",
                entity_id=contrast.cohort_id,
                path="excluded_genes",
                message=(
                    f"Excluded {len(excluded_gene_audits)} of {total_input_genes} input gene rows from canonical DE analysis due to harmonization constraints."
                ),
            )
        )

    valid_entries.sort(key=lambda e: e[0])
    valid_canonicals = [e[0] for e in valid_entries]
    valid_orig_ids = [e[1] for e in valid_entries]
    valid_row_indices = [e[2] for e in valid_entries]

    if not valid_canonicals:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="no_eligible_canonical_genes",
                entity_type="gene_harmonization",
                entity_id=contrast.cohort_id,
                path="canonical_genes",
                message="No uniquely mapped canonical genes found in raw count matrix.",
            )
        )
        return _make_blocked_result(
            contrast,
            findings,
            config_version,
            de_version,
            "No eligible canonical genes.",
            excluded_genes=tuple(excluded_gene_audits),
            total_input_genes_count=total_input_genes,
        )

    # 6. Design Resolution via DesignResolver
    design_res = DesignResolver.resolve(contrast, sample_metadata)
    findings.extend(design_res.findings)

    if design_res.status != ContrastStatus.ELIGIBLE or not design_res.is_estimable:
        return DifferentialExpressionContrastResult(
            contrast_id=contrast.contrast_id,
            study_id=contrast.study_id,
            experiment_id=contrast.experiment_id,
            cohort_id=contrast.cohort_id,
            relationship_id=contrast.relationship_id,
            treatment_condition_id=contrast.treatment_condition_id,
            matched_control_condition_ids=contrast.matched_control_condition_ids,
            treatment_sample_ids=contrast.treatment_sample_ids,
            control_sample_ids=contrast.control_sample_ids,
            treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
            control_biological_replicate_ids=contrast.control_biological_replicate_ids,
            treatment_replicate_count=contrast.treatment_replicate_count,
            control_replicate_count=contrast.control_replicate_count,
            treatment_exposures=contrast.treatment_exposures,
            sample_level_exposure_deviations=contrast.sample_level_exposure_deviations,
            agent_name=contrast.agent_name,
            agent_identifier=contrast.agent_identifier,
            vehicle=contrast.vehicle,
            concentration_or_dose=contrast.concentration_or_dose,
            concentration_or_dose_unit=contrast.concentration_or_dose_unit,
            exposure_start_time_or_stage=contrast.exposure_start_time_or_stage,
            developmental_age_or_stage_at_exposure=contrast.developmental_age_or_stage_at_exposure,
            exposure_duration=contrast.exposure_duration,
            exposure_duration_unit=contrast.exposure_duration_unit,
            washout_or_recovery_duration=contrast.washout_or_recovery_duration,
            washout_or_recovery_duration_unit=contrast.washout_or_recovery_duration_unit,
            treatment_collection_age_or_stage=contrast.treatment_collection_age_or_stage,
            control_collection_age_or_stage=contrast.control_collection_age_or_stage,
            treatment_collection_age_normalized=contrast.treatment_collection_age_normalized,
            control_collection_age_normalized=contrast.control_collection_age_normalized,
            biological_source_ids=contrast.biological_source_ids,
            organoid_context_ids=contrast.organoid_context_ids,
            design_formula=design_res.formula,
            design_rationale=design_res.rationale,
            design_matrix_rank=design_res.rank,
            residual_degrees_of_freedom=design_res.residual_degrees_of_freedom,
            size_factors_used=size_factors_used,
            excluded_genes=tuple(excluded_gene_audits),
            total_input_genes_count=total_input_genes,
            status=design_res.status,
            findings=ordered_findings(findings),
            configuration_version=config_version,
            de_version=de_version,
        )

    # 7. Build derived DESeq2 input tables
    # Ordering: Controls first (sorted by sample_id), then Treatments (sorted by sample_id)
    # This guarantees deterministic sample ordering regardless of input ordering
    ordered_ctrl = sorted(contrast.control_sample_ids)
    ordered_trt = sorted(contrast.treatment_sample_ids)
    ordered_samples = tuple(ordered_ctrl + ordered_trt)

    # Resolve biological source mapping for coldata if needed
    source_map: dict[str, str] = {}
    if isinstance(sample_metadata, dict):
        for sid, meta in sample_metadata.items():
            if isinstance(meta, dict):
                src = meta.get("biological_source_id") or meta.get("source_id") or meta.get("cell_line_id")
                if src:
                    source_map[sid] = str(src).strip()
    elif isinstance(sample_metadata, Sequence):
        for item in sample_metadata:
            if isinstance(item, dict) and "sample_id" in item:
                src = item.get("biological_source_id") or item.get("source_id") or item.get("cell_line_id")
                if src:
                    source_map[item["sample_id"]] = str(src).strip()

    coldata_rows: list[dict[str, str]] = []
    for sid in ordered_samples:
        cond = "control" if sid in ordered_ctrl else "treatment"
        row: dict[str, str] = {"sample_id": sid, "condition": cond}
        if "biological_source" in design_res.formula:
            row["biological_source"] = source_map.get(sid, contrast.biological_source_ids[0] if contrast.biological_source_ids else "source_1")
        coldata_rows.append(row)

    # Extract raw count columns for contrast samples
    # Mapping sample_id -> column in raw_matrix
    raw_sample_idx_map = {sid: idx for idx, sid in enumerate(raw_matrix.sample_ids)}

    count_cols: list[tuple[int, ...]] = []
    for sid in ordered_samples:
        col_idx = raw_sample_idx_map[sid]
        raw_col = raw_matrix.columns[col_idx]
        subset_col = tuple(raw_col[r_idx] for r_idx in valid_row_indices)
        count_cols.append(subset_col)

    # 8. Call official DESeq2 execution bridge
    r_result, env_info, r_findings, r_prov = execute_deseq2_contrast_r(
        gene_ids=tuple(valid_canonicals),
        sample_ids=ordered_samples,
        count_columns=tuple(count_cols),
        coldata=coldata_rows,
        size_factors=size_factors_used,
        design_formula=design_res.formula,
        r_binary_path=r_binary_path,
        strict_version_check=strict_version_check,
        r_timeout_seconds=r_timeout_seconds,
    )
    findings.extend(r_findings)

    if r_result is None or any(f.severity == Severity.ERROR for f in r_findings):
        return DifferentialExpressionContrastResult(
            contrast_id=contrast.contrast_id,
            study_id=contrast.study_id,
            experiment_id=contrast.experiment_id,
            cohort_id=contrast.cohort_id,
            relationship_id=contrast.relationship_id,
            treatment_condition_id=contrast.treatment_condition_id,
            matched_control_condition_ids=contrast.matched_control_condition_ids,
            treatment_sample_ids=contrast.treatment_sample_ids,
            control_sample_ids=contrast.control_sample_ids,
            treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
            control_biological_replicate_ids=contrast.control_biological_replicate_ids,
            treatment_replicate_count=contrast.treatment_replicate_count,
            control_replicate_count=contrast.control_replicate_count,
            treatment_exposures=contrast.treatment_exposures,
            sample_level_exposure_deviations=contrast.sample_level_exposure_deviations,
            agent_name=contrast.agent_name,
            agent_identifier=contrast.agent_identifier,
            vehicle=contrast.vehicle,
            concentration_or_dose=contrast.concentration_or_dose,
            concentration_or_dose_unit=contrast.concentration_or_dose_unit,
            exposure_start_time_or_stage=contrast.exposure_start_time_or_stage,
            developmental_age_or_stage_at_exposure=contrast.developmental_age_or_stage_at_exposure,
            exposure_duration=contrast.exposure_duration,
            exposure_duration_unit=contrast.exposure_duration_unit,
            washout_or_recovery_duration=contrast.washout_or_recovery_duration,
            washout_or_recovery_duration_unit=contrast.washout_or_recovery_duration_unit,
            treatment_collection_age_or_stage=contrast.treatment_collection_age_or_stage,
            control_collection_age_or_stage=contrast.control_collection_age_or_stage,
            treatment_collection_age_normalized=contrast.treatment_collection_age_normalized,
            control_collection_age_normalized=contrast.control_collection_age_normalized,
            biological_source_ids=contrast.biological_source_ids,
            organoid_context_ids=contrast.organoid_context_ids,
            design_formula=design_res.formula,
            design_rationale=design_res.rationale,
            design_matrix_rank=design_res.rank,
            residual_degrees_of_freedom=design_res.residual_degrees_of_freedom,
            size_factors_used=size_factors_used,
            excluded_genes=tuple(excluded_gene_audits),
            total_input_genes_count=total_input_genes,
            status=ContrastStatus.BLOCKED,
            findings=ordered_findings(findings),
            r_environment_info=env_info.to_dict(),
            input_hashes=r_prov.get("input_hashes", {}),
            configuration_version=config_version,
            de_version=de_version,
        )

    # 9. Format gene-wise output models
    gene_results: list[GeneDifferentialExpressionResult] = []
    for g_dict in r_result.get("gene_results", []):
        cid = g_dict["canonical_gene_id"]
        bm = float(g_dict["base_mean"]) if g_dict.get("base_mean") is not None else 0.0
        lfc = float(g_dict["log2_fold_change"]) if g_dict.get("log2_fold_change") is not None else None
        se = float(g_dict["lfc_standard_error"]) if g_dict.get("lfc_standard_error") is not None else None
        st = float(g_dict["wald_statistic"]) if g_dict.get("wald_statistic") is not None else None
        pv = float(g_dict["p_value"]) if g_dict.get("p_value") is not None else None
        padj = float(g_dict["adjusted_p_value_bh"]) if g_dict.get("adjusted_p_value_bh") is not None else None
        stat_str = g_dict.get("result_status", "OK")
        notes_str = g_dict.get("result_notes")

        gene_results.append(
            GeneDifferentialExpressionResult(
                canonical_gene_id=cid,
                approved_symbol=symbol_map.get(cid),
                base_mean=bm,
                log2_fold_change=lfc,
                lfc_standard_error=se,
                wald_statistic=st,
                p_value=pv,
                adjusted_p_value_bh=padj,
                result_status=stat_str,
                result_notes=notes_str,
            )
        )

    # Deterministic sorting of gene results by canonical_gene_id
    gene_results.sort(key=lambda x: x.canonical_gene_id)

    return DifferentialExpressionContrastResult(
        contrast_id=contrast.contrast_id,
        study_id=contrast.study_id,
        experiment_id=contrast.experiment_id,
        cohort_id=contrast.cohort_id,
        relationship_id=contrast.relationship_id,
        treatment_condition_id=contrast.treatment_condition_id,
        matched_control_condition_ids=contrast.matched_control_condition_ids,
        treatment_sample_ids=contrast.treatment_sample_ids,
        control_sample_ids=contrast.control_sample_ids,
        treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
        control_biological_replicate_ids=contrast.control_biological_replicate_ids,
        treatment_replicate_count=contrast.treatment_replicate_count,
        control_replicate_count=contrast.control_replicate_count,
        treatment_exposures=contrast.treatment_exposures,
        sample_level_exposure_deviations=contrast.sample_level_exposure_deviations,
        agent_name=contrast.agent_name,
        agent_identifier=contrast.agent_identifier,
        vehicle=contrast.vehicle,
        concentration_or_dose=contrast.concentration_or_dose,
        concentration_or_dose_unit=contrast.concentration_or_dose_unit,
        exposure_start_time_or_stage=contrast.exposure_start_time_or_stage,
        developmental_age_or_stage_at_exposure=contrast.developmental_age_or_stage_at_exposure,
        exposure_duration=contrast.exposure_duration,
        exposure_duration_unit=contrast.exposure_duration_unit,
        washout_or_recovery_duration=contrast.washout_or_recovery_duration,
        washout_or_recovery_duration_unit=contrast.washout_or_recovery_duration_unit,
        treatment_collection_age_or_stage=contrast.treatment_collection_age_or_stage,
        control_collection_age_or_stage=contrast.control_collection_age_or_stage,
        treatment_collection_age_normalized=contrast.treatment_collection_age_normalized,
        control_collection_age_normalized=contrast.control_collection_age_normalized,
        biological_source_ids=contrast.biological_source_ids,
        organoid_context_ids=contrast.organoid_context_ids,
        design_formula=design_res.formula,
        design_rationale=design_res.rationale,
        design_matrix_rank=r_result.get("design_matrix_rank", design_res.rank),
        residual_degrees_of_freedom=r_result.get("residual_degrees_of_freedom", design_res.residual_degrees_of_freedom),
        size_factors_used=size_factors_used,
        gene_results=tuple(gene_results),
        excluded_genes=tuple(excluded_gene_audits),
        total_input_genes_count=total_input_genes,
        status=ContrastStatus.ELIGIBLE,
        findings=ordered_findings(findings),
        r_environment_info=env_info.to_dict(),
        input_hashes=r_prov.get("input_hashes", {}),
        configuration_version=config_version,
        de_version=de_version,
    )


def run_differential_expression(
    contrast_dataset: ResponseContrastDataset,
    normalized_dataset: NormalizedDataset,
    harmonized_dataset: HarmonizedDataset,
    sample_metadata: dict[str, Any] | Sequence[dict[str, Any]] | None = None,
    r_binary_path: str = "Rscript",
    strict_version_check: bool = True,
    r_timeout_seconds: int = 180,
    config_version: str = "1.0.0",
    de_version: str = "1.0.0",
) -> DifferentialExpressionDataset:
    """Execute differential expression across all contrasts in a ResponseContrastDataset."""
    dataset_findings: list[Finding] = list(contrast_dataset.findings)
    results: list[DifferentialExpressionContrastResult] = []

    for contrast in contrast_dataset.contrasts:
        res = run_contrast_differential_expression(
            contrast=contrast,
            normalized_dataset=normalized_dataset,
            harmonized_dataset=harmonized_dataset,
            sample_metadata=sample_metadata,
            r_binary_path=r_binary_path,
            strict_version_check=strict_version_check,
            r_timeout_seconds=r_timeout_seconds,
            config_version=config_version,
            de_version=de_version,
        )
        results.append(res)
        dataset_findings.extend(res.findings)

    agg_status = status_for(dataset_findings)

    return DifferentialExpressionDataset(
        status=agg_status,
        findings=ordered_findings(dataset_findings),
        contrast_results=tuple(results),
        source_asset_id=contrast_dataset.source_asset_id,
        reference_identity=contrast_dataset.reference_identity,
        builder_version=contrast_dataset.builder_version,
        de_version=de_version,
        execution_provenance={
            "number_of_contrasts": len(results),
            "number_of_eligible": sum(1 for r in results if r.status == ContrastStatus.ELIGIBLE),
            "number_of_blocked": sum(1 for r in results if r.status == ContrastStatus.BLOCKED),
            "number_of_needs_review": sum(1 for r in results if r.status == ContrastStatus.NEEDS_REVIEW),
        },
    )


def _make_blocked_result(
    contrast: MolecularResponseContrast,
    findings: list[Finding],
    config_version: str,
    de_version: str,
    rationale: str,
    excluded_genes: tuple[ExcludedGeneAudit, ...] = (),
    total_input_genes_count: int = 0,
) -> DifferentialExpressionContrastResult:
    return DifferentialExpressionContrastResult(
        contrast_id=contrast.contrast_id,
        study_id=contrast.study_id,
        experiment_id=contrast.experiment_id,
        cohort_id=contrast.cohort_id,
        relationship_id=contrast.relationship_id,
        treatment_condition_id=contrast.treatment_condition_id,
        matched_control_condition_ids=contrast.matched_control_condition_ids,
        treatment_sample_ids=contrast.treatment_sample_ids,
        control_sample_ids=contrast.control_sample_ids,
        treatment_biological_replicate_ids=contrast.treatment_biological_replicate_ids,
        control_biological_replicate_ids=contrast.control_biological_replicate_ids,
        treatment_replicate_count=contrast.treatment_replicate_count,
        control_replicate_count=contrast.control_replicate_count,
        treatment_exposures=contrast.treatment_exposures,
        sample_level_exposure_deviations=contrast.sample_level_exposure_deviations,
        agent_name=contrast.agent_name,
        agent_identifier=contrast.agent_identifier,
        vehicle=contrast.vehicle,
        concentration_or_dose=contrast.concentration_or_dose,
        concentration_or_dose_unit=contrast.concentration_or_dose_unit,
        exposure_start_time_or_stage=contrast.exposure_start_time_or_stage,
        developmental_age_or_stage_at_exposure=contrast.developmental_age_or_stage_at_exposure,
        exposure_duration=contrast.exposure_duration,
        exposure_duration_unit=contrast.exposure_duration_unit,
        washout_or_recovery_duration=contrast.washout_or_recovery_duration,
        washout_or_recovery_duration_unit=contrast.washout_or_recovery_duration_unit,
        treatment_collection_age_or_stage=contrast.treatment_collection_age_or_stage,
        control_collection_age_or_stage=contrast.control_collection_age_or_stage,
        treatment_collection_age_normalized=contrast.treatment_collection_age_normalized,
        control_collection_age_normalized=contrast.control_collection_age_normalized,
        biological_source_ids=contrast.biological_source_ids,
        organoid_context_ids=contrast.organoid_context_ids,
        design_formula="",
        design_rationale=rationale,
        excluded_genes=excluded_genes,
        total_input_genes_count=total_input_genes_count,
        status=ContrastStatus.BLOCKED,
        findings=ordered_findings(findings),
        configuration_version=config_version,
        de_version=de_version,
    )
