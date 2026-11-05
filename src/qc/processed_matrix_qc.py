"""Quality control for processed continuous bulk RNA-seq expression matrices.

Implements Phase U3.2 Step 6C Processed Matrix QC:
1. Validates continuous floating-point expression matrices without imposing raw-integer-count requirements.
2. Evaluates continuous sample-level metrics (total expression, detected genes, zero fraction, quantiles, min, median, max).
3. Reconciles matrix sample headers against canonical metadata assays and samples.
4. Preserves duplicate gene rows and un-aggregated biological replication.
5. Evaluates within-condition biological replicate rank correlations deterministically.
6. Enforces shared QC severity and status precedence models.
"""

from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from itertools import combinations
from math import sqrt
from pathlib import Path
from statistics import median
from typing import Any

from src.ingestion.processed_matrix_reader import ProcessedExpressionMatrix, read_processed_matrix
from src.ingestion.registry import reported_value
from src.qc.metrics import quantile
from src.qc.result import Finding, PairCorrelation, Severity, Status, ordered, status_for

PROCESSED_QC_RULESET_VERSION = "processed_matrix_qc_v1"
PROCESSED_QC_CALCULATION_VERSION = "processed_matrix_qc_engine_v1"


def ranks(values: Sequence[float]) -> list[float]:
    """Compute standard average ranks for ties deterministically."""
    out = [0.0] * len(values)
    ordered_items = sorted(enumerate(values), key=lambda x: x[1])
    i = 0
    n = len(ordered_items)
    while i < n:
        j = i
        while j + 1 < n and ordered_items[j + 1][1] == ordered_items[i][1]:
            j += 1
        rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            out[ordered_items[k][0]] = rank
        i = j + 1
    return out


def spearman_rank_correlation(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Compute Spearman rank correlation between two continuous numeric vectors."""
    if len(a) != len(b) or len(a) < 2:
        return None
    ra = ranks(a)
    rb = ranks(b)
    ma = sum(ra) / len(ra)
    mb = sum(rb) / len(rb)
    da = [x - ma for x in ra]
    db = [x - mb for x in rb]
    denominator = sqrt(sum(x * x for x in da) * sum(x * x for x in db))
    if denominator == 0.0:
        return None
    return sum(x * y for x, y in zip(da, db)) / denominator


@dataclass(frozen=True, slots=True)
class ProcessedSampleMetrics:
    """Descriptive statistics for a single sample in a processed expression matrix."""

    sample_id: str
    condition_id: str | None
    biological_replicate_id: str | None
    technical_replicate_id: str | None
    total_expression: float
    detected_genes: int
    zero_expression_genes: int
    zero_fraction: float
    minimum: float
    median: float
    maximum: float
    mean: float
    quantiles: tuple[float, float, float]
    findings: tuple[Finding, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "condition_id": self.condition_id,
            "biological_replicate_id": self.biological_replicate_id,
            "technical_replicate_id": self.technical_replicate_id,
            "total_expression": self.total_expression,
            "detected_genes": self.detected_genes,
            "zero_expression_genes": self.zero_expression_genes,
            "zero_fraction": self.zero_fraction,
            "minimum": self.minimum,
            "median": self.median,
            "maximum": self.maximum,
            "mean": self.mean,
            "quantiles": list(self.quantiles),
            "findings": [asdict(f) for f in self.findings],
        }


@dataclass(frozen=True, slots=True)
class ProcessedConditionMetrics:
    """Condition-level summary of biological replicate rank correlations for continuous matrices."""

    condition_id: str
    sample_count: int
    biological_replicate_count: int
    technical_replicates: tuple[tuple[str, tuple[str, ...]], ...]
    correlation_method: str
    correlation_implementation: str
    pairs: tuple[PairCorrelation, ...]
    evaluable_pairs: int
    correlation_summary: tuple[float, float, float] | None
    unevaluable_reasons: tuple[str, ...]
    findings: tuple[Finding, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "condition_id": self.condition_id,
            "sample_count": self.sample_count,
            "biological_replicate_count": self.biological_replicate_count,
            "technical_replicates": [list(item) for item in self.technical_replicates],
            "correlation_method": self.correlation_method,
            "correlation_implementation": self.correlation_implementation,
            "pairs": [asdict(p) for p in self.pairs],
            "evaluable_pairs": self.evaluable_pairs,
            "correlation_summary": list(self.correlation_summary) if self.correlation_summary else None,
            "unevaluable_reasons": list(self.unevaluable_reasons),
            "findings": [asdict(f) for f in self.findings],
        }


@dataclass(frozen=True, slots=True)
class ProcessedMatrixQCResult:
    """Deterministic quality control outcome for a continuous processed expression matrix."""

    status: Status
    findings: tuple[Finding, ...]
    number_of_genes: int
    number_of_samples: int
    source_asset_id: str
    assay_id: str
    study_ids: tuple[str, ...]
    experiment_ids: tuple[str, ...]
    matrix: ProcessedExpressionMatrix | None
    samples: tuple[ProcessedSampleMetrics, ...]
    conditions: tuple[ProcessedConditionMetrics, ...]
    ruleset_version: str
    calculation_version: str

    @property
    def is_valid(self) -> bool:
        """True if the QC evaluation passed without fatal errors."""
        return self.status != Status.FAIL

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "findings": [asdict(f) for f in self.findings],
            "number_of_genes": self.number_of_genes,
            "number_of_samples": self.number_of_samples,
            "source_asset_id": self.source_asset_id,
            "assay_id": self.assay_id,
            "study_ids": list(self.study_ids),
            "experiment_ids": list(self.experiment_ids),
            "samples": [s.to_dict() for s in self.samples],
            "conditions": [c.to_dict() for c in self.conditions],
            "ruleset_version": self.ruleset_version,
            "calculation_version": self.calculation_version,
        }


def reconcile_processed_metadata(
    metadata: dict[str, Any],
    matrix_samples: Sequence[str],
    asset_id: str | None,
    assay_id: str | None,
) -> tuple[dict[str, dict[str, Any]], list[Finding]]:
    """Reconcile sample headers against canonical metadata without inferring biological identities."""
    findings: list[Finding] = []
    samples = [s for s in metadata.get("samples", []) if isinstance(s, dict)]
    index: dict[str, list[dict[str, Any]]] = {}
    for s in samples:
        sid = s.get("sample_id")
        if sid:
            index.setdefault(sid, []).append(s)

    assay = None
    if assay_id:
        assay = next((a for a in metadata.get("sequencing_assays", []) if a.get("assay_id") == assay_id), None)
    elif metadata.get("sequencing_assays"):
        assay = metadata.get("sequencing_assays")[0]
        assay_id = assay.get("assay_id")

    asset = None
    if asset_id:
        asset = next((a for a in metadata.get("input_data_assets", []) if a.get("asset_id") == asset_id), None)
    elif metadata.get("input_data_assets"):
        asset = metadata.get("input_data_assets")[0]
        asset_id = asset.get("asset_id")

    expected: set[str] = set()
    if assay:
        expected = set(assay.get("sample_ids", []))
    elif assay_id:
        findings.append(Finding(Severity.ERROR, "assay_missing", "assay", assay_id, "assay_id", "Assay not found in canonical metadata."))

    if asset and assay_id:
        asset_assays = asset.get("assay_ids", [])
        if assay_id not in asset_assays:
            findings.append(Finding(Severity.ERROR, "asset_assay_mismatch", "asset", asset_id, "assay_ids", f"Asset '{asset_id}' is not linked to assay '{assay_id}'."))

    for sid in matrix_samples:
        if sid not in index:
            findings.append(Finding(Severity.REVIEW, "unexpected_matrix_samples", "sample", sid, "matrix.header", f"Matrix sample '{sid}' is absent from canonical metadata."))
        elif len(index[sid]) != 1:
            findings.append(Finding(Severity.ERROR, "ambiguous_sample_mapping", "sample", sid, "samples", f"Sample '{sid}' does not map uniquely in canonical metadata."))
        elif expected and sid not in expected:
            findings.append(Finding(Severity.ERROR, "assay_sample_mismatch", "sample", sid, "assay.sample_ids", f"Matrix sample '{sid}' is not assigned to assay '{assay_id}'."))

    for sid in sorted(expected - set(matrix_samples)):
        findings.append(Finding(Severity.REVIEW, "missing_expected_samples", "sample", sid, "matrix.header", f"Expected assay sample '{sid}' is absent from expression matrix."))

    contexts: dict[str, dict[str, Any]] = {}
    for sid in matrix_samples:
        matching = index.get(sid)
        if matching and len(matching) == 1:
            s = matching[0]
            contexts[sid] = {
                "condition_id": s.get("condition_id"),
                "experiment_id": s.get("experiment_id"),
                "biological_replicate_id": reported_value(s.get("biological_replicate_id")),
                "technical_replicate_id": reported_value(s.get("technical_replicate_id")),
            }
        else:
            contexts[sid] = {
                "condition_id": None,
                "experiment_id": None,
                "biological_replicate_id": None,
                "technical_replicate_id": None,
            }

    return contexts, findings


def validate_processed_matrix(
    matrix_input: ProcessedExpressionMatrix | Path | str,
    metadata: dict[str, Any],
    asset_id: str | None = None,
    assay_id: str | None = None,
    config: dict[str, Any] | None = None,
) -> ProcessedMatrixQCResult:
    """Perform deterministic QC validation of a processed continuous expression matrix.

    Evaluates structural integrity, numerical validity, sample reconciliation,
    and condition-level biological replicate correlation without imposing raw-count integer constraints.
    """
    all_findings: list[Finding] = []

    # 1. Resolve matrix input
    matrix: ProcessedExpressionMatrix | None = None
    if isinstance(matrix_input, (str, Path)):
        intake_data = read_processed_matrix(matrix_input)
        for f in intake_data.findings:
            sev = Severity[f.severity] if f.severity in Severity.__members__ else Severity.ERROR
            all_findings.append(Finding(sev, f.rule, "matrix", None, f.path, f.message))
        if not intake_data.is_valid or intake_data.matrix is None:
            # Fatal structural intake failure
            ordered_findings = ordered(all_findings)
            return ProcessedMatrixQCResult(
                status=Status.FAIL,
                findings=ordered_findings,
                number_of_genes=0,
                number_of_samples=0,
                source_asset_id=asset_id or "",
                assay_id=assay_id or "",
                study_ids=tuple(s.get("study_id", "") for s in metadata.get("studies", [])),
                experiment_ids=tuple(e.get("experiment_id", "") for e in metadata.get("experiments", [])),
                matrix=None,
                samples=(),
                conditions=(),
                ruleset_version=PROCESSED_QC_RULESET_VERSION,
                calculation_version=PROCESSED_QC_CALCULATION_VERSION,
            )
        matrix = intake_data.matrix
    else:
        matrix = matrix_input

    # 2. Extract context and study/experiment IDs
    studies = tuple(s.get("study_id", "") for s in metadata.get("studies", []))
    experiments = tuple(e.get("experiment_id", "") for e in metadata.get("experiments", []))
    resolved_asset_id = asset_id or (metadata.get("input_data_assets", [{}])[0].get("asset_id", "") if metadata.get("input_data_assets") else "")
    resolved_assay_id = assay_id or (metadata.get("sequencing_assays", [{}])[0].get("assay_id", "") if metadata.get("sequencing_assays") else "")

    # 3. Canonical metadata reconciliation
    contexts, recon_findings = reconcile_processed_metadata(
        metadata=metadata,
        matrix_samples=matrix.sample_ids,
        asset_id=resolved_asset_id,
        assay_id=resolved_assay_id,
    )
    all_findings.extend(recon_findings)

    # 4. Compute continuous sample-level metrics
    sample_results: list[ProcessedSampleMetrics] = []
    for s_idx, sample_id in enumerate(matrix.sample_ids):
        col = matrix.columns[s_idx]
        ctx = contexts.get(sample_id, {})
        s_findings: list[Finding] = []

        total_expr = sum(col)
        zero_genes = sum(1 for v in col if v == 0.0)
        detected = len(col) - zero_genes
        zero_fraction = zero_genes / len(col) if col else 0.0
        min_val = min(col) if col else 0.0
        median_val = median(col) if col else 0.0
        max_val = max(col) if col else 0.0
        mean_val = total_expr / len(col) if col else 0.0
        quantiles = (
            quantile(col, 0.25),
            quantile(col, 0.50),
            quantile(col, 0.75),
        ) if col else (0.0, 0.0, 0.0)

        if col and all(v == 0.0 for v in col):
            s_finding = Finding(
                Severity.REVIEW,
                "zero_only_samples",
                "sample",
                sample_id,
                f"sample[{sample_id}]",
                f"Sample '{sample_id}' has zero expression across all genes.",
            )
            s_findings.append(s_finding)
            all_findings.append(s_finding)

        sample_results.append(
            ProcessedSampleMetrics(
                sample_id=sample_id,
                condition_id=ctx.get("condition_id"),
                biological_replicate_id=ctx.get("biological_replicate_id"),
                technical_replicate_id=ctx.get("technical_replicate_id"),
                total_expression=total_expr,
                detected_genes=detected,
                zero_expression_genes=zero_genes,
                zero_fraction=zero_fraction,
                minimum=min_val,
                median=median_val,
                maximum=max_val,
                mean=mean_val,
                quantiles=quantiles,
                findings=tuple(s_findings),
            )
        )

    # 5. Dataset-level gene observations
    for g_idx, gene_id in enumerate(matrix.gene_ids):
        if matrix.number_of_samples > 0 and all(matrix.columns[s_idx][g_idx] == 0.0 for s_idx in range(matrix.number_of_samples)):
            all_findings.append(
                Finding(
                    Severity.INFO,
                    "zero_only_genes",
                    "gene",
                    gene_id,
                    f"gene[{gene_id}]",
                    f"Gene '{gene_id}' has zero expression in every sample.",
                )
            )

    # 6. Condition-level biological replicate correlation
    condition_groups: dict[str | None, list[tuple[ProcessedSampleMetrics, tuple[float, ...]]]] = {}
    for res, s_idx in zip(sample_results, range(matrix.number_of_samples)):
        condition_groups.setdefault(res.condition_id, []).append((res, matrix.columns[s_idx]))

    condition_results: list[ProcessedConditionMetrics] = []
    for cond_id, items in sorted(condition_groups.items(), key=lambda x: x[0] or ""):
        if not cond_id:
            continue
        by_rep: dict[str | None, list[tuple[ProcessedSampleMetrics, tuple[float, ...]]]] = {}
        technical: dict[str, list[str]] = {}

        for res, col in items:
            by_rep.setdefault(res.biological_replicate_id, []).append((res, col))
            if res.technical_replicate_id is not None and res.biological_replicate_id is not None:
                technical.setdefault(res.biological_replicate_id, []).append(res.technical_replicate_id)

        pairs: list[PairCorrelation] = []
        reasons: list[str] = []
        cond_findings: list[Finding] = []

        rep_keys = sorted(by_rep.keys(), key=lambda x: x or "")
        for (ra, aa), (rb, bb) in combinations([(k, by_rep[k]) for k in rep_keys], 2):
            if len(aa) != 1 or len(bb) != 1:
                val = None
                reason = "technical_replicates_require_unapproved_aggregation"
                reasons.append(reason)
                f_item = Finding(
                    Severity.REVIEW,
                    reason,
                    "condition",
                    cond_id,
                    "technical_replicates",
                    "Technical replicates present without approved aggregation policy.",
                )
                cond_findings.append(f_item)
                all_findings.append(f_item)
            else:
                val = spearman_rank_correlation(aa[0][1], bb[0][1])
                reason = None if val is not None else "constant_or_zero_variance_vector"
                if reason:
                    reasons.append(reason)
            pairs.append(
                PairCorrelation(
                    sample_a=aa[0][0].sample_id,
                    sample_b=bb[0][0].sample_id,
                    replicate_a=ra or "",
                    replicate_b=rb or "",
                    value=val,
                    reason=reason,
                )
            )

        evaluable = [p.value for p in pairs if p.value is not None]
        summary = (min(evaluable), median(evaluable), max(evaluable)) if evaluable else None
        if len(by_rep) < 2:
            reasons.append("fewer_than_two_biological_replicates")

        condition_results.append(
            ProcessedConditionMetrics(
                condition_id=cond_id,
                sample_count=len(items),
                biological_replicate_count=len(by_rep),
                technical_replicates=tuple((k, tuple(sorted(v))) for k, v in sorted(technical.items())),
                correlation_method="spearman",
                correlation_implementation="spearman_rank_correlation",
                pairs=tuple(pairs),
                evaluable_pairs=len(evaluable),
                correlation_summary=summary,
                unevaluable_reasons=tuple(sorted(set(reasons))),
                findings=tuple(cond_findings),
            )
        )

    ordered_findings_tuple = ordered(all_findings)
    overall_status = status_for(ordered_findings_tuple)

    return ProcessedMatrixQCResult(
        status=overall_status,
        findings=ordered_findings_tuple,
        number_of_genes=matrix.number_of_genes,
        number_of_samples=matrix.number_of_samples,
        source_asset_id=resolved_asset_id,
        assay_id=resolved_assay_id,
        study_ids=studies,
        experiment_ids=experiments,
        matrix=matrix,
        samples=tuple(sample_results),
        conditions=tuple(condition_results),
        ruleset_version=PROCESSED_QC_RULESET_VERSION,
        calculation_version=PROCESSED_QC_CALCULATION_VERSION,
    )
