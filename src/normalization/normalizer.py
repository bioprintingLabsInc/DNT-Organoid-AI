"""Orchestration of Bulk RNA-seq Normalization v1 for human bulk expression data."""

from collections import defaultdict
import math
from pathlib import Path
from typing import Any

from src.gene_harmonization.models import (
    GeneMappingStatus,
    HarmonizedDataset,
    QCReviewDisposition,
    Status as HarmonizationStatus,
)
from src.qc.matrix import CountMatrix
from src.qc.result import Status as QCStatus

from .cohort import NormalizationCohortBuilder
from .config import NormalizationConfig, config_checksum, load_config, validate_config
from .models import (
    CohortNormalizationMetrics,
    CollisionDisposition,
    Finding,
    NormalizationCohort,
    NormalizationMethod,
    NormalizationReviewDisposition,
    NormalizedCohortResult,
    NormalizedCountMatrix,
    NormalizedDataset,
    ResponseEligibility,
    SampleSizeFactor,
    Severity,
    Status,
    ordered,
    status_for,
)
from .r_bridge import REnvironmentInfo, estimate_size_factors_r

NORMALIZATION_VERSION = "bulk_normalization_v1"


class BulkNormalizer:
    """Orchestrates metadata-derived cohort formation, official DESeq2 size-factor estimation,

    diagnostic normalized count calculation, and normalization technical quality control.
    """

    def __init__(self, config: NormalizationConfig | None = None) -> None:
        self.config = config or NormalizationConfig()

    def normalize(
        self,
        harmonized: HarmonizedDataset,
        metadata: dict[str, Any],
        config: NormalizationConfig | None = None,
        review_disposition: NormalizationReviewDisposition | dict[str, Any] | str | None = None,
        collision_disposition: CollisionDisposition | dict[str, Any] | str | None = None,
    ) -> NormalizedDataset:
        """Execute Bulk RNA-seq Normalization v1 across all cohorts represented in the dataset."""
        cfg = config or self.config
        validate_config(cfg.to_dict())
        findings: list[Finding] = []

        rev_disp = NormalizationReviewDisposition.from_any(review_disposition)
        col_disp = CollisionDisposition.from_any(collision_disposition)

        raw_matrix = harmonized.matrix
        source_asset_id = harmonized.source_asset_id

        # 1. Upstream QC Gating
        upstream_qc = harmonized.upstream_qc_status
        if upstream_qc == QCStatus.FAIL:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="upstream_qc_failure",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="upstream_qc_status",
                    message="Upstream raw-count QC failed. Normalization is completely blocked.",
                )
            )
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        if upstream_qc == QCStatus.NEEDS_REVIEW:
            if not harmonized.qc_review_disposition:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="unauthorized_upstream_qc_needs_review",
                        entity_type="dataset",
                        entity_id=source_asset_id,
                        path="qc_review_disposition",
                        message=(
                            "Upstream QC resulted in NEEDS_REVIEW but lacks an auditable QCReviewDisposition. "
                            "Normalization is blocked."
                        ),
                    )
                )
                return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)
            findings.append(
                Finding(
                    severity=Severity.INFO,
                    rule_id="upstream_qc_needs_review_authorized",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="qc_review_disposition",
                    message=f"Proceeding under authorized upstream QC review disposition: {harmonized.qc_review_disposition.decision}.",
                )
            )

        # 2. Upstream Harmonization Gating
        if harmonized.status == HarmonizationStatus.FAIL:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="upstream_harmonization_failure",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="harmonization_status",
                    message="Upstream Gene Harmonization failed. Normalization is blocked.",
                )
            )
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        if raw_matrix is None:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="missing_raw_count_matrix",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="matrix",
                    message="Input HarmonizedDataset has no count matrix.",
                )
            )
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        # 3. Canonical Collision Gate
        if harmonized.collisions:
            if not col_disp:
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="canonical_collisions_detected",
                        entity_type="dataset",
                        entity_id=source_asset_id,
                        path="collisions",
                        message=(
                            f"Input dataset contains {len(harmonized.collisions)} canonical collision groups. "
                            "Default behavior blocks normalization at NEEDS_REVIEW with no automatic collapsing. "
                            "Requires an explicit auditable CollisionDisposition."
                        ),
                    )
                )
                # Raw matrix is preserved; normalized output is unavailable (None)
                return self._build_empty_dataset(
                    Status.NEEDS_REVIEW,
                    findings,
                    source_asset_id,
                    harmonized,
                    cfg,
                    collision_disposition=col_disp,
                )
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="canonical_collisions_authorized",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="collisions",
                    message=f"Proceeding past collision gate under authorized CollisionDisposition: {col_disp.decision}.",
                )
            )

        # 4. Feature Eligibility
        # Filter strictly for UNIQUELY_MAPPED rows with non-null canonical_gene_id and no collisions
        colliding_indices: set[int] = set()
        for col in harmonized.collisions:
            colliding_indices.update(col.source_indices)

        eligible_indices: list[int] = []
        for i, g in enumerate(harmonized.genes):
            if g.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED and g.canonical_gene_id and i not in colliding_indices:
                eligible_indices.append(i)

        if not eligible_indices:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="no_eligible_features",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="genes",
                    message="No features in the harmonized dataset are eligible for size-factor estimation.",
                )
            )
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        # 5. Cohort Formation
        cohort_builder = NormalizationCohortBuilder(metadata)
        cohorts, cohort_findings = cohort_builder.build_cohorts(harmonized.sample_ids, source_asset_id)
        findings.extend(cohort_findings)

        if any(f.severity == Severity.ERROR for f in cohort_findings):
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        if not cohorts:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="no_cohorts_constructed",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="cohorts",
                    message="Failed to construct any valid Normalization Cohort from matrix samples.",
                )
            )
            return self._build_empty_dataset(Status.FAIL, findings, source_asset_id, harmonized, cfg)

        # 6. Normalize Each Cohort
        cohort_results: list[NormalizedCohortResult] = []
        sample_to_col_idx = {sid: idx for idx, sid in enumerate(raw_matrix.sample_ids)}

        for cohort in cohorts:
            cohort_res = self._normalize_single_cohort(
                cohort=cohort,
                raw_matrix=raw_matrix,
                eligible_indices=eligible_indices,
                harmonized_genes=harmonized.genes,
                cohort_builder=cohort_builder,
                config=cfg,
                review_disposition=rev_disp,
                collision_disposition=col_disp,
            )
            cohort_results.append(cohort_res)
            findings.extend(cohort_res.findings)

        overall_status = status_for(findings)
        return NormalizedDataset(
            status=overall_status,
            findings=ordered(findings),
            cohort_results=tuple(cohort_results),
            source_asset_id=source_asset_id,
            reference_identity=harmonized.reference_identity,
            config_version=cfg.config_version,
            normalization_version=NORMALIZATION_VERSION,
        )

    def _normalize_single_cohort(
        self,
        cohort: NormalizationCohort,
        raw_matrix: CountMatrix,
        eligible_indices: list[int],
        harmonized_genes: tuple[Any, ...],
        cohort_builder: NormalizationCohortBuilder,
        config: NormalizationConfig,
        review_disposition: NormalizationReviewDisposition | None,
        collision_disposition: CollisionDisposition | None,
    ) -> NormalizedCohortResult:
        findings: list[Finding] = []
        method = config.size_factor_method

        # Method authorization check for poscounts
        if method == NormalizationMethod.POSCOUNTS:
            is_authorized = (
                review_disposition is not None and "poscounts" in review_disposition.decision.lower()
            ) or (config.size_factor_method == "poscounts")
            if not is_authorized:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="unauthorized_poscounts_requested",
                        entity_type="cohort",
                        entity_id=cohort.cohort_id,
                        path="size_factor_method",
                        message=(
                            "type='poscounts' requested without explicit scientific authorization. "
                            "Automatic or unreviewed execution of poscounts is strictly prohibited."
                        ),
                    )
                )
                return self._build_empty_cohort_result(
                    Status.FAIL, cohort, raw_matrix, findings, method, review_disposition, collision_disposition
                )

        # Extract eligible counts submatrix for cohort samples
        sample_indices = [raw_matrix.sample_ids.index(sid) for sid in cohort.sample_ids]
        eligible_columns: list[tuple[int, ...]] = []
        for s_idx in sample_indices:
            full_col = raw_matrix.columns[s_idx]
            elig_col = tuple(full_col[i] for i in eligible_indices)
            eligible_columns.append(elig_col)

        eligible_gene_ids = tuple(harmonized_genes[i].canonical_gene_id for i in eligible_indices)

        # Invoke Official DESeq2 R execution bridge
        sfs_dict, env_info, r_findings, provenance = estimate_size_factors_r(
            matrix_columns=tuple(eligible_columns),
            gene_ids=eligible_gene_ids,
            sample_ids=cohort.sample_ids,
            method=method,
            config=config,
        )
        findings.extend(r_findings)

        if sfs_dict is None:
            # If standard ratio zero-geometric-mean failure:
            is_standard_zero_fail = any(f.rule_id == "standard_size_factors_failed" for f in r_findings)
            cohort_status = Status.NEEDS_REVIEW if is_standard_zero_fail else Status.FAIL

            return self._build_empty_cohort_result(
                cohort_status,
                cohort,
                raw_matrix,
                findings,
                method,
                review_disposition,
                collision_disposition,
                env_info=env_info.to_dict(),
                provenance=provenance,
            )

        # Size factors successfully estimated!
        # Resolve matched controls and ResponseEligibility
        sample_controls = cohort_builder.resolve_sample_controls(cohort)

        # Compute raw library sizes and normalized count columns
        total_genes = len(raw_matrix.gene_ids)
        norm_columns: list[tuple[float, ...]] = []
        sample_size_factor_objs: list[SampleSizeFactor] = []

        for sid in cohort.sample_ids:
            s_idx = raw_matrix.sample_ids.index(sid)
            raw_col = raw_matrix.columns[s_idx]
            raw_lib_size = sum(raw_col)
            sf = sfs_dict[sid]

            # Scale ALL rows in raw matrix by sf for diagnostic/inspection matrix
            norm_col = tuple(float(val) / sf for val in raw_col)
            norm_columns.append(norm_col)
            norm_lib_size = sum(norm_col)

            elig, matched_controls, tc_status = sample_controls.get(
                sid, (ResponseEligibility.UNASSIGNED_OR_AMBIGUOUS, (), "unassigned")
            )

            sample_size_factor_objs.append(
                SampleSizeFactor(
                    sample_id=sid,
                    size_factor=sf,
                    raw_library_size=raw_lib_size,
                    normalized_library_size=norm_lib_size,
                    cohort_id=cohort.cohort_id,
                    response_eligibility=elig,
                    matched_control_sample_ids=matched_controls,
                    treatment_control_status=tc_status,
                )
            )

        diagnostic_norm_matrix = NormalizedCountMatrix(
            gene_ids=raw_matrix.gene_ids,
            sample_ids=cohort.sample_ids,
            columns=tuple(norm_columns),
        )

        # Calculate descriptive QC metrics (without invented thresholds)
        all_sfs = [sf_obj.size_factor for sf_obj in sample_size_factor_objs]
        sf_min = min(all_sfs)
        sf_max = max(all_sfs)
        sf_median = float(sorted(all_sfs)[len(all_sfs) // 2])
        sf_mean = float(sum(all_sfs) / len(all_sfs))
        sf_variance = float(sum((x - sf_mean) ** 2 for x in all_sfs) / len(all_sfs)) if len(all_sfs) > 1 else 0.0
        sf_std = math.sqrt(sf_variance)
        sf_ratio = sf_max / sf_min if sf_min > 0 else float("inf")

        # Count positive geometric mean features across cohort samples
        m = len(cohort.sample_ids)
        pos_geo_count = 0
        total_cells = total_genes * m
        zero_cells = 0

        for row_idx in eligible_indices:
            # Gene is positive geometric mean in standard ratio if counts > 0 in all cohort samples
            has_zero = any(raw_matrix.columns[s_idx][row_idx] == 0 for s_idx in sample_indices)
            if not has_zero:
                pos_geo_count += 1

        for s_idx in sample_indices:
            zero_cells += raw_matrix.columns[s_idx].count(0)

        zero_fraction = float(zero_cells / total_cells) if total_cells > 0 else 0.0

        n_treat = len(cohort.treatment_sample_ids)
        n_ctrl = len(cohort.control_sample_ids)
        tc_ratio = float(n_treat / n_ctrl) if n_ctrl > 0 else None

        metrics = CohortNormalizationMetrics(
            total_features_count=total_genes,
            eligible_features_count=len(eligible_indices),
            positive_geometric_mean_features_count=pos_geo_count,
            percentage_features_evaluable=float(pos_geo_count / total_genes * 100.0) if total_genes > 0 else 0.0,
            size_factor_min=sf_min,
            size_factor_max=sf_max,
            size_factor_median=sf_median,
            size_factor_mean=sf_mean,
            size_factor_std=sf_std,
            size_factor_ratio_max_min=sf_ratio,
            number_of_samples=m,
            number_of_treatment_samples=n_treat,
            number_of_control_samples=n_ctrl,
            treatment_control_ratio=tc_ratio,
            zero_count_fraction=zero_fraction,
        )

        # Check configurable thresholds (if any non-null approved threshold is configured)
        if config.min_evaluable_genes is not None and metrics.eligible_features_count < config.min_evaluable_genes:
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="configured_threshold_exceeded",
                    entity_type="cohort",
                    entity_id=cohort.cohort_id,
                    path="metrics.eligible_features_count",
                    message=(
                        f"Eligible features count ({metrics.eligible_features_count}) is below configured threshold "
                        f"({config.min_evaluable_genes})."
                    ),
                )
            )

        if config.max_size_factor_ratio is not None and sf_ratio > config.max_size_factor_ratio:
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="configured_threshold_exceeded",
                    entity_type="cohort",
                    entity_id=cohort.cohort_id,
                    path="metrics.size_factor_ratio_max_min",
                    message=(
                        f"Max/min size factor ratio ({sf_ratio:.2f}) exceeds configured threshold "
                        f"({config.max_size_factor_ratio:.2f})."
                    ),
                )
            )

        cohort_status = status_for(findings)

        return NormalizedCohortResult(
            status=cohort_status,
            cohort=cohort,
            findings=ordered(findings),
            size_factors=tuple(sample_size_factor_objs),
            metrics=metrics,
            diagnostic_normalized_matrix=diagnostic_norm_matrix,
            raw_count_matrix=raw_matrix,
            method=method,
            r_environment_info=env_info.to_dict(),
            execution_provenance=provenance,
            review_disposition=review_disposition,
            collision_disposition=collision_disposition,
        )

    def _build_empty_dataset(
        self,
        status: Status,
        findings: list[Finding],
        source_asset_id: str,
        harmonized: HarmonizedDataset,
        config: NormalizationConfig,
        collision_disposition: CollisionDisposition | None = None,
    ) -> NormalizedDataset:
        return NormalizedDataset(
            status=status,
            findings=ordered(findings),
            cohort_results=(),
            source_asset_id=source_asset_id,
            reference_identity=harmonized.reference_identity,
            config_version=config.config_version,
            normalization_version=NORMALIZATION_VERSION,
        )

    def _build_empty_cohort_result(
        self,
        status: Status,
        cohort: NormalizationCohort,
        raw_matrix: CountMatrix,
        findings: list[Finding],
        method: str,
        review_disposition: NormalizationReviewDisposition | None,
        collision_disposition: CollisionDisposition | None,
        env_info: dict[str, Any] | None = None,
        provenance: dict[str, Any] | None = None,
    ) -> NormalizedCohortResult:
        return NormalizedCohortResult(
            status=status,
            cohort=cohort,
            findings=ordered(findings),
            size_factors=(),
            metrics=None,
            diagnostic_normalized_matrix=None,
            raw_count_matrix=raw_matrix,
            method=method,
            r_environment_info=env_info or {},
            execution_provenance=provenance or {},
            review_disposition=review_disposition,
            collision_disposition=collision_disposition,
        )


def normalize_bulk_dataset(
    harmonized: HarmonizedDataset,
    metadata: dict[str, Any],
    config: NormalizationConfig | None = None,
    review_disposition: NormalizationReviewDisposition | dict[str, Any] | str | None = None,
    collision_disposition: CollisionDisposition | dict[str, Any] | str | None = None,
) -> NormalizedDataset:
    """Convenience function executing Bulk RNA-seq Normalization v1."""
    normalizer = BulkNormalizer(config=config)
    return normalizer.normalize(
        harmonized=harmonized,
        metadata=metadata,
        config=config,
        review_disposition=review_disposition,
        collision_disposition=collision_disposition,
    )
