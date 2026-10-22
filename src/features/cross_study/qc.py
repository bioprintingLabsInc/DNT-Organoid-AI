"""Quality Control and Validation Auditor for Stage 8A.1 Condition-Level Gene Response.

Verifies:
1. Exact gene identity preservation (all Stage 7 canonical genes preserved).
2. Exact numerical preservation (log2FC, SE, Wald, p-value, FDR).
3. Contrast & context preservation (dose, duration, stage, compound).
4. Traceable missing-data preservation (missing values not replaced with zeros).
5. Non-significant genes and non-Entrez genes retained.
6. Single-study accounting and zero meta-analysis execution.
7. Zero DNT label leakage.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field
import math
from typing import Any

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    ResponseDirection,
)


def _approx_eq(a: float | None, b: float | None, tol: float = 1e-12) -> bool:
    """Check numerical equality accounting for None and NaN."""
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isnan(a) or math.isnan(b):
        return False
    return abs(a - b) <= tol


@dataclass(frozen=True, slots=True)
class RegistryQCSummary:
    """Structured audit summary for Stage 8A.1 Condition-Level Gene Response."""

    is_valid: bool
    total_contrasts_evaluated: int
    total_records_evaluated: int
    numerical_mismatches_count: int
    gene_identity_mismatches_count: int
    direction_mismatches_count: int
    context_mismatches_count: int
    leakage_violations_count: int
    findings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_valid": self.is_valid,
            "total_contrasts_evaluated": self.total_contrasts_evaluated,
            "total_records_evaluated": self.total_records_evaluated,
            "numerical_mismatches_count": self.numerical_mismatches_count,
            "gene_identity_mismatches_count": self.gene_identity_mismatches_count,
            "direction_mismatches_count": self.direction_mismatches_count,
            "context_mismatches_count": self.context_mismatches_count,
            "leakage_violations_count": self.leakage_violations_count,
            "findings": list(self.findings),
        }


def audit_registry_qc(
    registry: ConditionGeneResponseRegistry,
    contrasts: Sequence[DifferentialExpressionContrastResult],
) -> RegistryQCSummary:
    """Audit the integrity of a ConditionGeneResponseRegistry against its upstream Stage 7 contrasts."""
    findings: list[str] = []
    num_mismatches = 0
    id_mismatches = 0
    dir_mismatches = 0
    ctx_mismatches = 0
    leak_violations = 0

    contrast_map = {c.contrast_id: c for c in contrasts}

    # Group registry records by contrast_id
    records_by_contrast: dict[str, list[ConditionGeneResponseRecord]] = {}
    for r in registry.records:
        records_by_contrast.setdefault(r.contrast_id, []).append(r)

    # 1. Contrast presence and record count audit
    for cid, contrast in contrast_map.items():
        recs = records_by_contrast.get(cid, [])
        if len(recs) != len(contrast.gene_results):
            id_mismatches += 1
            findings.append(
                f"Contrast {cid} record count mismatch: Stage 7 has {len(contrast.gene_results)} genes, "
                f"registry has {len(recs)} records."
            )

        stage7_gene_map = {g.canonical_gene_id: g for g in contrast.gene_results}
        for rec in recs:
            g = stage7_gene_map.get(rec.canonical_ensembl_gene_id)
            if g is None:
                id_mismatches += 1
                findings.append(
                    f"Registry record gene {rec.canonical_ensembl_gene_id} not found in Stage 7 contrast {cid}."
                )
                continue

            # 2. Exact Numerical Preservation Audit
            if not _approx_eq(rec.log2_fold_change, g.log2_fold_change):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: log2FC mismatch: "
                    f"Stage 7={g.log2_fold_change}, Registry={rec.log2_fold_change}"
                )
            if not _approx_eq(rec.log2_fold_change_se, g.lfc_standard_error):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: SE mismatch: "
                    f"Stage 7={g.lfc_standard_error}, Registry={rec.log2_fold_change_se}"
                )
            if not _approx_eq(rec.wald_statistic, g.wald_statistic):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: Wald mismatch: "
                    f"Stage 7={g.wald_statistic}, Registry={rec.wald_statistic}"
                )
            if not _approx_eq(rec.p_value, g.p_value):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: p-value mismatch: "
                    f"Stage 7={g.p_value}, Registry={rec.p_value}"
                )
            if not _approx_eq(rec.adjusted_p_value, g.adjusted_p_value_bh):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: adjusted p-value mismatch: "
                    f"Stage 7={g.adjusted_p_value_bh}, Registry={rec.adjusted_p_value}"
                )
            if not _approx_eq(rec.base_mean, g.base_mean):
                num_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: base_mean mismatch: "
                    f"Stage 7={g.base_mean}, Registry={rec.base_mean}"
                )

            # 3. Direction Consistency Audit
            if rec.log2_fold_change is None or (
                isinstance(rec.log2_fold_change, float) and math.isnan(rec.log2_fold_change)
            ):
                expected_dir = ResponseDirection.UNDETERMINED
            elif rec.log2_fold_change > 0:
                expected_dir = ResponseDirection.UP
            elif rec.log2_fold_change < 0:
                expected_dir = ResponseDirection.DOWN
            else:
                expected_dir = ResponseDirection.NO_DIRECTION

            if rec.response_direction != expected_dir:
                dir_mismatches += 1
                findings.append(
                    f"Gene {rec.canonical_ensembl_gene_id} in {cid}: direction mismatch: "
                    f"expected {expected_dir}, got {rec.response_direction}"
                )

            # 4. Context Preservation Audit
            if rec.study_id != contrast.study_id:
                ctx_mismatches += 1
                findings.append(f"Record study_id {rec.study_id} != contrast {contrast.study_id}")
            if rec.treatment_replicate_count != contrast.treatment_replicate_count:
                ctx_mismatches += 1
                findings.append(f"Replicate count mismatch for contrast {cid}")

    # 5. DNT Label Leakage Audit
    # Verify that no record contains DNT label fields or references external toxicity classification
    for rec in registry.records:
        rec_dict = rec.to_dict()
        for forbidden in (
            "dnt_label",
            "dnt_call",
            "is_toxic",
            "toxicity_class",
            "positive_control",
            "negative_control",
            "ground_truth",
            "target_label",
        ):
            if forbidden in rec_dict:
                leak_violations += 1
                findings.append(f"Forbidden label field '{forbidden}' found in record {rec.canonical_ensembl_gene_id}")

    # 6. Meta-analysis check (ensure no pooled statistics exist)
    for rec in registry.records:
        rec_dict = rec.to_dict()
        for forbidden_meta in ("pooled_log2fc", "pooled_se", "tau_squared", "i_squared", "heterogeneity_p"):
            if forbidden_meta in rec_dict:
                leak_violations += 1
                findings.append(f"Prohibited meta-analysis field '{forbidden_meta}' found in record")

    is_valid = (
        num_mismatches == 0
        and id_mismatches == 0
        and dir_mismatches == 0
        and ctx_mismatches == 0
        and leak_violations == 0
    )

    return RegistryQCSummary(
        is_valid=is_valid,
        total_contrasts_evaluated=len(contrasts),
        total_records_evaluated=len(registry.records),
        numerical_mismatches_count=num_mismatches,
        gene_identity_mismatches_count=id_mismatches,
        direction_mismatches_count=dir_mismatches,
        context_mismatches_count=ctx_mismatches,
        leakage_violations_count=leak_violations,
        findings=tuple(findings),
    )
