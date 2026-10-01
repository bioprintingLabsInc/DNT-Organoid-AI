"""Deterministic DesignResolver for Step 7B Differential Expression v1.

Evaluates sample metadata, biological replication, biological source structure,
and estimability to choose the correct statistical model formula or block execution.
"""

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from src.response_builder.models import ContrastStatus, Finding, MolecularResponseContrast, Severity


def matrix_rank(matrix: Sequence[Sequence[float]], tol: float = 1e-9) -> int:
    """Compute matrix rank via Gaussian elimination with partial pivoting."""
    a = [list(map(float, row)) for row in matrix]
    rows = len(a)
    if rows == 0:
        return 0
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot_row = None
        max_val = tol
        for r in range(rank, rows):
            if abs(a[r][col]) > max_val:
                max_val = abs(a[r][col])
                pivot_row = r
        if pivot_row is None:
            continue
        a[rank], a[pivot_row] = a[pivot_row], a[rank]
        pivot = a[rank][col]
        for r in range(rows):
            if r != rank and abs(a[r][col]) > tol:
                factor = a[r][col] / pivot
                for c in range(col, cols):
                    a[r][c] -= factor * a[rank][c]
        rank += 1
    return rank


@dataclass(frozen=True, slots=True)
class DesignResolutionResult:
    """Outcome of model design resolution and estimability validation."""

    formula: str
    rationale: str
    status: ContrastStatus
    findings: tuple[Finding, ...]
    rank: int | None = None
    residual_degrees_of_freedom: int | None = None
    is_estimable: bool = False
    samples_used: tuple[str, ...] = ()


class DesignResolver:
    """Scientific design resolver for differential expression models."""

    @classmethod
    def resolve(
        cls,
        contrast: MolecularResponseContrast,
        sample_metadata: dict[str, Any] | Sequence[dict[str, Any]] | None = None,
    ) -> DesignResolutionResult:
        """Resolve the model design formula and evaluate estimability.

        Rules:
        - Case A: Single biological source / cell line with independent biological replicates -> ~ condition
        - Case B: Multiple paired biological sources represented in both arms -> ~ biological_source + condition
        - Case C: Independent unpaired unique biological sources (each replicate has unique source ID) -> ~ condition
        - Case D: Treatment is inseparable from biological source -> BLOCKED (treatment_source_confounding)
        - Case E: Partially crossed / ambiguous biological sources -> NEEDS_REVIEW (partially_crossed_biological_sources)
        - Replication: Minimum 2 distinct biological replicates in treatment and control arms
        - Technical replicates: Multiple sample columns with identical biological_replicate_id trigger NEEDS_REVIEW
        - Estimability: Model matrix must be full rank, treatment effect estimable, residual df > 0
        """
        findings: list[Finding] = []

        trt_samples = tuple(sorted(contrast.treatment_sample_ids))
        ctrl_samples = tuple(sorted(contrast.control_sample_ids))
        all_samples = ctrl_samples + trt_samples  # control first as reference level

        if not trt_samples or not ctrl_samples:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="missing_contrast_samples",
                    entity_type="contrast",
                    entity_id=contrast.contrast_id,
                    path="samples",
                    message="Contrast must contain both treatment and control samples.",
                )
            )
            return DesignResolutionResult(
                formula="",
                rationale="Missing contrast samples.",
                status=ContrastStatus.BLOCKED,
                findings=tuple(findings),
            )

        # Build sample lookup map
        meta_by_id: dict[str, dict[str, Any]] = {}
        if isinstance(sample_metadata, dict):
            for k, v in sample_metadata.items():
                if isinstance(v, dict):
                    meta_by_id[k] = v
                else:
                    meta_by_id[k] = {"sample_id": k}
        elif isinstance(sample_metadata, Sequence):
            for item in sample_metadata:
                if isinstance(item, dict) and "sample_id" in item:
                    meta_by_id[item["sample_id"]] = item

        # Helper to get bio replicate ID
        def get_bio_rep(sid: str, is_trt: bool) -> str:
            if sid in meta_by_id and meta_by_id[sid].get("biological_replicate_id"):
                return str(meta_by_id[sid]["biological_replicate_id"]).strip()
            # Fallback to contrast bio replicate tuples if 1:1 length
            if is_trt:
                idx = contrast.treatment_sample_ids.index(sid) if sid in contrast.treatment_sample_ids else -1
                if 0 <= idx < len(contrast.treatment_biological_replicate_ids):
                    return str(contrast.treatment_biological_replicate_ids[idx]).strip()
            else:
                idx = contrast.control_sample_ids.index(sid) if sid in contrast.control_sample_ids else -1
                if 0 <= idx < len(contrast.control_biological_replicate_ids):
                    return str(contrast.control_biological_replicate_ids[idx]).strip()
            return sid

        # Helper to get bio source ID
        def get_bio_source(sid: str) -> str | None:
            if sid in meta_by_id:
                m = meta_by_id[sid]
                src = (
                    m.get("biological_source_id")
                    or m.get("source_id")
                    or m.get("cell_line_id")
                    or m.get("donor_id")
                )
                if src:
                    return str(src).strip()
            if len(contrast.biological_source_ids) == 1:
                return contrast.biological_source_ids[0]
            return None

        # 1. Check technical replicates and biological replication
        trt_bio_reps = [get_bio_rep(s, is_trt=True) for s in trt_samples]
        ctrl_bio_reps = [get_bio_rep(s, is_trt=False) for s in ctrl_samples]

        trt_rep_counts = Counter(trt_bio_reps)
        ctrl_rep_counts = Counter(ctrl_bio_reps)

        has_technical_reps = any(c > 1 for c in trt_rep_counts.values()) or any(c > 1 for c in ctrl_rep_counts.values())
        if has_technical_reps:
            dup_reps = [r for r, c in trt_rep_counts.items() if c > 1] + [r for r, c in ctrl_rep_counts.items() if c > 1]
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="unresolved_technical_replicates",
                    entity_type="contrast",
                    entity_id=contrast.contrast_id,
                    path="biological_replicate_ids",
                    message=(
                        f"Multiple count-matrix sample columns share biological_replicate_id(s) {dup_reps} "
                        "without an approved technical-replicate aggregation rule."
                    ),
                )
            )
            return DesignResolutionResult(
                formula="",
                rationale="Unresolved technical replicates detected; not passed independently to DESeq2.",
                status=ContrastStatus.NEEDS_REVIEW,
                findings=tuple(findings),
            )

        unique_trt_reps = set(trt_bio_reps)
        unique_ctrl_reps = set(ctrl_bio_reps)

        if len(unique_trt_reps) < 2 or len(unique_ctrl_reps) < 2:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="insufficient_biological_replicates",
                    entity_type="contrast",
                    entity_id=contrast.contrast_id,
                    path="biological_replicate_count",
                    message=(
                        f"Insufficient biological replicates: treatment has {len(unique_trt_reps)}, "
                        f"control has {len(unique_ctrl_reps)} (minimum 2 required in each arm for differential expression inference)."
                    ),
                )
            )
            return DesignResolutionResult(
                formula="",
                rationale="Insufficient biological replicates in treatment or control arm.",
                status=ContrastStatus.BLOCKED,
                findings=tuple(findings),
            )

        # 2. Analyze biological sources
        trt_sources = {get_bio_source(s) for s in trt_samples}
        ctrl_sources = {get_bio_source(s) for s in ctrl_samples}
        all_sources_set = trt_sources | ctrl_sources

        # Clean None values if present
        clean_trt_sources = {s for s in trt_sources if s is not None}
        clean_ctrl_sources = {s for s in ctrl_sources if s is not None}
        clean_all_sources = sorted(clean_trt_sources | clean_ctrl_sources)

        formula = "~ condition"
        rationale = ""
        include_source_covariate = False

        if len(clean_all_sources) <= 1:
            # Case A: Single biological source / cell line
            formula = "~ condition"
            rationale = "Single biological source; condition is the only varying factor."
        else:
            # Multiple sources present
            is_disjoint = clean_trt_sources.isdisjoint(clean_ctrl_sources)
            if is_disjoint:
                # Check whether source IDs repeat within condition
                all_source_assignments = [get_bio_source(s) for s in all_samples]
                source_counts = Counter([s for s in all_source_assignments if s is not None])
                has_repeating_source = any(cnt > 1 for cnt in source_counts.values())

                if has_repeating_source:
                    # Case D: Confounded treatment and biological source
                    findings.append(
                        Finding(
                            severity=Severity.ERROR,
                            rule_id="treatment_source_confounding",
                            entity_type="contrast",
                            entity_id=contrast.contrast_id,
                            path="biological_source_ids",
                            message=(
                                f"Treatment condition is completely confounded with biological source: "
                                f"treatment sources {sorted(clean_trt_sources)} vs control sources {sorted(clean_ctrl_sources)}."
                            ),
                        )
                    )
                    return DesignResolutionResult(
                        formula="",
                        rationale="Treatment condition is completely confounded with biological source.",
                        status=ContrastStatus.BLOCKED,
                        findings=tuple(findings),
                    )
                else:
                    # Case C: Independent unpaired unique biological sources (each sample has unique source ID)
                    formula = "~ condition"
                    rationale = (
                        "Independent unpaired biological sources where each biological source represents an "
                        "independent biological replicate without reuse; analyzed using ~ condition."
                    )
            elif clean_trt_sources == clean_ctrl_sources:
                # Case B: Multiple paired / blocked biological sources represented in both arms
                # Verify each source has both treatment and control
                paired_all = all(
                    any(get_bio_source(s) == src for s in trt_samples)
                    and any(get_bio_source(s) == src for s in ctrl_samples)
                    for src in clean_all_sources
                )
                if paired_all:
                    formula = "~ biological_source + condition"
                    rationale = "Multiple paired biological sources represented in both arms; biological source included as blocking factor."
                    include_source_covariate = True
                else:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="partially_crossed_biological_sources",
                            entity_type="contrast",
                            entity_id=contrast.contrast_id,
                            path="biological_source_ids",
                            message=(
                                f"Biological sources are partially crossed between treatment ({sorted(clean_trt_sources)}) "
                                f"and control ({sorted(clean_ctrl_sources)}) arms; requires expert review."
                            ),
                        )
                    )
                    return DesignResolutionResult(
                        formula="",
                        rationale="Partially crossed biological source structure requires expert review.",
                        status=ContrastStatus.NEEDS_REVIEW,
                        findings=tuple(findings),
                    )
            else:
                # Case E: Partially crossed / ambiguous biological sources
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="partially_crossed_biological_sources",
                        entity_type="contrast",
                        entity_id=contrast.contrast_id,
                        path="biological_source_ids",
                        message=(
                            f"Partially crossed biological sources: treatment ({sorted(clean_trt_sources)}) "
                            f"vs control ({sorted(clean_ctrl_sources)}) arms require expert review."
                        ),
                    )
                )
                return DesignResolutionResult(
                    formula="",
                    rationale="Partially crossed biological sources require expert review.",
                    status=ContrastStatus.NEEDS_REVIEW,
                    findings=tuple(findings),
                )

        # 3. Construct model matrix and evaluate rank, degrees of freedom, estimability
        x_rows: list[list[float]] = []
        n_samples = len(all_samples)

        if include_source_covariate:
            # Columns: [intercept, (S-1 dummy columns for sources 1..S-1), condition_treatment]
            ref_source = clean_all_sources[0]
            other_sources = clean_all_sources[1:]
            for s in all_samples:
                row = [1.0]
                s_source = get_bio_source(s)
                for src in other_sources:
                    row.append(1.0 if s_source == src else 0.0)
                is_trt = 1.0 if s in trt_samples else 0.0
                row.append(is_trt)
                x_rows.append(row)
        else:
            # Columns: [intercept, condition_treatment]
            for s in all_samples:
                is_trt = 1.0 if s in trt_samples else 0.0
                x_rows.append([1.0, is_trt])

        p_cols = len(x_rows[0])
        rank = matrix_rank(x_rows)
        resid_df = n_samples - rank

        # Check rank deficiency
        if rank < p_cols:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="rank_deficient_design",
                    entity_type="model_design",
                    entity_id=formula,
                    path="design_matrix",
                    message=f"Model matrix is rank deficient: rank {rank} < {p_cols} columns.",
                )
            )
            return DesignResolutionResult(
                formula=formula,
                rationale=rationale,
                status=ContrastStatus.BLOCKED,
                findings=tuple(findings),
                rank=rank,
                residual_degrees_of_freedom=resid_df,
                is_estimable=False,
            )

        # Check condition estimability
        x_no_cond = [row[:-1] for row in x_rows]
        rank_no_cond = matrix_rank(x_no_cond)
        if rank_no_cond == rank:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="treatment_effect_not_estimable",
                    entity_type="model_design",
                    entity_id=formula,
                    path="condition_effect",
                    message="Treatment condition is completely collinear with other design covariates; treatment effect is not estimable.",
                )
            )
            return DesignResolutionResult(
                formula=formula,
                rationale=rationale,
                status=ContrastStatus.BLOCKED,
                findings=tuple(findings),
                rank=rank,
                residual_degrees_of_freedom=resid_df,
                is_estimable=False,
            )

        # Check residual degrees of freedom
        if resid_df <= 0:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="insufficient_residual_degrees_of_freedom",
                    entity_type="model_design",
                    entity_id=formula,
                    path="residual_degrees_of_freedom",
                    message=f"Residual degrees of freedom must be > 0 (got {resid_df} for {n_samples} samples and rank {rank}).",
                )
            )
            return DesignResolutionResult(
                formula=formula,
                rationale=rationale,
                status=ContrastStatus.BLOCKED,
                findings=tuple(findings),
                rank=rank,
                residual_degrees_of_freedom=resid_df,
                is_estimable=False,
            )

        return DesignResolutionResult(
            formula=formula,
            rationale=rationale,
            status=ContrastStatus.ELIGIBLE,
            findings=tuple(findings),
            rank=rank,
            residual_degrees_of_freedom=resid_df,
            is_estimable=True,
            samples_used=all_samples,
        )
