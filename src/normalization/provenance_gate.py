"""Processing provenance gate and scientific compatibility validator for processed expression matrices.

Implements Phase U3.2 Step 6C Provenance Gate:
1. Evaluates documented measurement units (TPM, FPKM, CPM, logCPM, normalized_counts).
2. Distinguishes normalization methods, transformation scales, and processing histories.
3. Strictly forbids guessing whether data is linear or log-transformed.
4. Emits INSUFFICIENT_PROCESSING_PROVENANCE and blocks inferential processing if essential provenance is absent.
5. Verifies physical and mathematical scale compatibility (e.g. non-negativity of linear or offset log scales).
6. Scientific isolation: strictly blocks DESeq2 routing and prevents rounding to artificial integer counts.
7. Produces validated U3.2 output container (ProcessedExpressionDataset) ready for future U3.3 limma-trend engine.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import StrEnum
import math
from pathlib import Path
from typing import Any

from src.gene_harmonization.models import GeneIdentifierType, GeneMappingStatus, HarmonizedGene
from src.gene_harmonization.reference import GeneReference
from src.ingestion.processed_matrix_reader import ProcessedExpressionMatrix, ProcessedMatrixData, read_processed_matrix
from src.qc.processed_matrix_qc import ProcessedMatrixQCResult, validate_processed_matrix
from src.qc.result import Finding, Severity, Status, ordered, status_for


class ExpressionUnit(StrEnum):
    """Recognized biological measurement / expression units."""

    TPM = "TPM"
    FPKM = "FPKM"
    CPM = "CPM"
    LOG_CPM = "logCPM"
    NORMALIZED_COUNTS = "normalized_counts"
    RAW_COUNTS = "raw_counts"
    UNKNOWN = "unknown"


class TransformationScale(StrEnum):
    """Recognized numerical transformation scales."""

    LINEAR = "linear"
    LOG2 = "log2"
    LOG10 = "log10"
    LN = "ln"
    UNKNOWN = "unknown"


class ProcessedNormalizationMethod(StrEnum):
    """Recognized library-depth / gene-length normalization methodologies."""

    TPM = "tpm"
    FPKM = "fpkm"
    CPM = "cpm"
    DESEQ2_SIZE_FACTORS = "deseq2_size_factors"
    EDGER_TMM = "edger_tmm"
    QUANTILE = "quantile"
    NONE = "none"
    OTHER = "other"
    UNKNOWN = "unknown"


NormalizationMethod = ProcessedNormalizationMethod


@dataclass(frozen=True, slots=True)
class ProcessingProvenance:
    """Documented upstream processing provenance for a bulk RNA-seq expression matrix."""

    expression_unit: ExpressionUnit
    transformation_scale: TransformationScale
    normalization_method: NormalizationMethod
    pseudocount: float | None = None
    is_log_transformed: bool = False
    software_pipeline: str | None = None
    pipeline_version: str | None = None
    reference_annotation: str | None = None
    genome_build: str | None = None
    notes: str | None = None
    raw_metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "expression_unit": str(self.expression_unit),
            "transformation_scale": str(self.transformation_scale),
            "normalization_method": str(self.normalization_method),
            "pseudocount": self.pseudocount,
            "is_log_transformed": self.is_log_transformed,
            "software_pipeline": self.software_pipeline,
            "pipeline_version": self.pipeline_version,
            "reference_annotation": self.reference_annotation,
            "genome_build": self.genome_build,
            "notes": self.notes,
        }


@dataclass(frozen=True, slots=True)
class ProvenanceGateResult:
    """Outcome of processing provenance evaluation and downstream compatibility gating."""

    status: Status
    findings: tuple[Finding, ...]
    provenance: ProcessingProvenance
    is_inferentially_eligible: bool
    deseq2_eligible: bool
    deseq2_ineligibility_reason: str
    limma_trend_eligible: bool
    target_analysis_engine: str | None
    blocked_analysis_engines: tuple[str, ...] = ("DESeq2",)
    evaluated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    @property
    def is_valid(self) -> bool:
        return self.status != Status.FAIL

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "findings": [asdict(f) for f in self.findings],
            "provenance": self.provenance.to_dict(),
            "is_inferentially_eligible": self.is_inferentially_eligible,
            "deseq2_eligible": self.deseq2_eligible,
            "deseq2_ineligibility_reason": self.deseq2_ineligibility_reason,
            "limma_trend_eligible": self.limma_trend_eligible,
            "target_analysis_engine": self.target_analysis_engine,
            "blocked_analysis_engines": list(self.blocked_analysis_engines),
            "evaluated_at": self.evaluated_at,
        }


@dataclass(frozen=True, slots=True)
class ProcessedExpressionDataset:
    """Canonical validated U3.2 dataset ready for consumption by future U3.3 analysis engine."""

    dataset_id: str
    matrix: ProcessedExpressionMatrix
    provenance: ProcessingProvenance
    qc_result: ProcessedMatrixQCResult
    gate_result: ProvenanceGateResult
    harmonized_genes: tuple[HarmonizedGene, ...]
    eligible_canonical_genes_count: int
    excluded_genes_count: int
    total_input_genes_count: int
    status: Status
    is_inferentially_eligible: bool
    target_analysis_engine: str | None
    blocked_analysis_engines: tuple[str, ...] = ("DESeq2",)

    @property
    def is_gene_universe_conserved(self) -> bool:
        """Law of Conservation of the Gene Universe: total == eligible + excluded."""
        return self.total_input_genes_count == (self.eligible_canonical_genes_count + self.excluded_genes_count)

    def to_dict(self) -> dict[str, Any]:
        return {
            "dataset_id": self.dataset_id,
            "status": str(self.status),
            "is_inferentially_eligible": self.is_inferentially_eligible,
            "target_analysis_engine": self.target_analysis_engine,
            "blocked_analysis_engines": list(self.blocked_analysis_engines),
            "total_input_genes_count": self.total_input_genes_count,
            "eligible_canonical_genes_count": self.eligible_canonical_genes_count,
            "excluded_genes_count": self.excluded_genes_count,
            "is_gene_universe_conserved": self.is_gene_universe_conserved,
            "provenance": self.provenance.to_dict(),
            "qc_status": str(self.qc_result.status),
            "gate_status": str(self.gate_result.status),
            "matrix_dimensions": {
                "number_of_genes": self.matrix.number_of_genes,
                "number_of_samples": self.matrix.number_of_samples,
            },
        }


def parse_processing_provenance(
    metadata: dict[str, Any],
    overrides: dict[str, Any] | None = None,
    asset_id: str | None = None,
) -> ProcessingProvenance:
    """Extract and standardize processing provenance from metadata without guessing."""
    raw = {}
    if asset_id and metadata.get("input_data_assets"):
        matched_asset = next((a for a in metadata["input_data_assets"] if a.get("asset_id") == asset_id), None)
        if matched_asset and isinstance(matched_asset.get("processing_provenance"), dict):
            raw.update(matched_asset["processing_provenance"])
        elif matched_asset:
            raw.update(matched_asset)

    if isinstance(metadata.get("processing_provenance"), dict):
        raw.update(metadata["processing_provenance"])

    if overrides:
        raw.update(overrides)

    # 1. Parse Expression Unit
    unit_str = str(
        raw.get("expression_unit")
        or raw.get("measurement_unit")
        or raw.get("unit")
        or ""
    ).strip().lower()

    if unit_str in {"tpm", "transcripts per million", "transcript per million"}:
        unit = ExpressionUnit.TPM
    elif unit_str in {"fpkm", "rpkm", "fragments per kilobase million"}:
        unit = ExpressionUnit.FPKM
    elif unit_str in {"cpm", "counts per million"}:
        unit = ExpressionUnit.CPM
    elif unit_str in {"logcpm", "log_cpm", "log2cpm", "log2_cpm", "log2-cpm"}:
        unit = ExpressionUnit.LOG_CPM
    elif unit_str in {
        "normalized_counts", "normalized counts", "size_factor_normalized",
        "deseq2_normalized", "norm_counts", "normalized",
    }:
        unit = ExpressionUnit.NORMALIZED_COUNTS
    elif unit_str in {"raw_counts", "raw counts", "counts", "unnormalized_counts"}:
        unit = ExpressionUnit.RAW_COUNTS
    elif not unit_str:
        unit = ExpressionUnit.UNKNOWN
    else:
        unit = ExpressionUnit.UNKNOWN

    # 2. Parse Transformation Scale
    scale_str = str(
        raw.get("transformation_scale")
        or raw.get("transformation")
        or raw.get("scale")
        or ""
    ).strip().lower()

    pseudocount_val = raw.get("pseudocount")
    if pseudocount_val is not None:
        try:
            pseudocount_val = float(pseudocount_val)
        except (ValueError, TypeError):
            pseudocount_val = None

    if unit == ExpressionUnit.LOG_CPM and not scale_str:
        # logCPM explicitly declares log scale
        scale = TransformationScale.LOG2
        is_log = True
    elif scale_str in {"linear", "untransformed", "raw_scale", "none", "natural"}:
        scale = TransformationScale.LINEAR
        is_log = False
    elif scale_str in {
        "log2", "log2(x+1)", "log2_plus_one", "log2p1", "log2(tpm+1)", "log2(fpkm+1)", "log2(cpm+1)",
    }:
        scale = TransformationScale.LOG2
        is_log = True
        if pseudocount_val is None and any(tag in scale_str for tag in {"+1", "plus_one", "p1"}):
            pseudocount_val = 1.0
    elif scale_str in {"log10", "log10(x+1)"}:
        scale = TransformationScale.LOG10
        is_log = True
        if pseudocount_val is None and "+1" in scale_str:
            pseudocount_val = 1.0
    elif scale_str in {"ln", "natural_log"}:
        scale = TransformationScale.LN
        is_log = True
    elif not scale_str:
        scale = TransformationScale.UNKNOWN
        is_log = False
    else:
        scale = TransformationScale.UNKNOWN
        is_log = False

    # 3. Parse Normalization Method
    norm_str = str(
        raw.get("normalization_method")
        or raw.get("normalization")
        or ""
    ).strip().lower()

    if norm_str in {"tpm", "transcripts per million"}:
        norm = NormalizationMethod.TPM
    elif norm_str in {"fpkm", "rpkm"}:
        norm = NormalizationMethod.FPKM
    elif norm_str in {"cpm", "counts per million"}:
        norm = NormalizationMethod.CPM
    elif norm_str in {"deseq2", "deseq2_size_factors", "median_of_ratios", "median of ratios", "poscounts"}:
        norm = NormalizationMethod.DESEQ2_SIZE_FACTORS
    elif norm_str in {"edger", "edger_tmm", "tmm"}:
        norm = NormalizationMethod.EDGER_TMM
    elif norm_str in {"quantile", "quantile_normalization"}:
        norm = NormalizationMethod.QUANTILE
    elif norm_str in {"none", "unnormalized"}:
        norm = NormalizationMethod.NONE
    elif norm_str:
        norm = NormalizationMethod.OTHER
    elif unit in {ExpressionUnit.TPM, ExpressionUnit.FPKM, ExpressionUnit.CPM, ExpressionUnit.LOG_CPM}:
        # In the absence of an explicit normalization method string, TPM/FPKM/CPM units define their normalization
        norm = NormalizationMethod[unit.value.replace("logCPM", "CPM")]
    else:
        norm = NormalizationMethod.UNKNOWN

    return ProcessingProvenance(
        expression_unit=unit,
        transformation_scale=scale,
        normalization_method=norm,
        pseudocount=pseudocount_val,
        is_log_transformed=is_log,
        software_pipeline=raw.get("software_pipeline") or raw.get("software") or raw.get("pipeline"),
        pipeline_version=raw.get("pipeline_version") or raw.get("software_version"),
        reference_annotation=raw.get("reference_annotation") or raw.get("annotation"),
        genome_build=raw.get("genome_build"),
        notes=raw.get("notes"),
        raw_metadata=raw,
    )


def evaluate_provenance_gate(
    matrix: ProcessedExpressionMatrix | None,
    metadata: dict[str, Any],
    asset_id: str | None = None,
    provenance_overrides: dict[str, Any] | None = None,
    qc_result: ProcessedMatrixQCResult | None = None,
) -> ProvenanceGateResult:
    """Evaluate processing provenance and determine scientific suitability for downstream analysis.

    Enforces:
    1. Rejection of unknown or undocumented expression units (INSUFFICIENT_PROCESSING_PROVENANCE).
    2. Zero guessing of transformation scale: linear vs log must be explicitly documented.
    3. Rejection of negative values in linear expression scales.
    4. Rejection of negative values where incompatible with declared pseudocount offsets.
    5. Complete DESeq2 isolation: processed matrices are never eligible for DESeq2 raw-count pipeline.
    6. Qualification for future U3.3 limma-trend engine.
    """
    findings: list[Finding] = []
    provenance = parse_processing_provenance(metadata, overrides=provenance_overrides, asset_id=asset_id)

    # 1. Essential Expression Unit Provenance
    if provenance.expression_unit == ExpressionUnit.UNKNOWN:
        findings.append(
            Finding(
                Severity.ERROR,
                "INSUFFICIENT_PROCESSING_PROVENANCE",
                "processing_provenance",
                None,
                "expression_unit",
                "Documented expression unit (e.g. TPM, FPKM, CPM, normalized_counts) is required and cannot be established.",
            )
        )
    elif provenance.expression_unit == ExpressionUnit.RAW_COUNTS:
        findings.append(
            Finding(
                Severity.ERROR,
                "INCOMPATIBLE_RAW_COUNTS",
                "processing_provenance",
                None,
                "expression_unit",
                "Raw integer counts detected in processed expression intake; raw counts must be routed through BULK_RAW_COUNT_PIPELINE.",
            )
        )

    # 2. Essential Transformation Scale Provenance
    if provenance.transformation_scale == TransformationScale.UNKNOWN:
        findings.append(
            Finding(
                Severity.ERROR,
                "INSUFFICIENT_PROCESSING_PROVENANCE",
                "processing_provenance",
                None,
                "transformation_scale",
                "Transformation scale (e.g. linear, log2) is undocumented. System will not guess transformation scale.",
            )
        )

    # 3. Normalization Methodology / Processing History
    if (
        provenance.normalization_method == NormalizationMethod.UNKNOWN
        and not provenance.software_pipeline
        and provenance.expression_unit == ExpressionUnit.NORMALIZED_COUNTS
    ):
        findings.append(
            Finding(
                Severity.ERROR,
                "INSUFFICIENT_PROCESSING_PROVENANCE",
                "processing_provenance",
                None,
                "normalization_method",
                "Normalization methodology (e.g. DESeq2 size factors, edgeR TMM) is undocumented for normalized count matrix.",
            )
        )

    # 4. Mathematical and Scale Value Checks
    if matrix is not None and matrix.number_of_genes > 0 and matrix.number_of_samples > 0:
        all_values = [v for col in matrix.columns for v in col]

        # Non-finite check
        if any(math.isnan(v) or math.isinf(v) for v in all_values):
            findings.append(
                Finding(
                    Severity.ERROR,
                    "NON_FINITE_EXPRESSION_VALUES",
                    "matrix",
                    None,
                    "values",
                    "Matrix contains non-finite expression values (NaN or Inf).",
                )
            )

        min_val = min(all_values)
        max_val = max(all_values)

        # Linear scale non-negativity check
        if provenance.transformation_scale == TransformationScale.LINEAR:
            if min_val < 0.0:
                findings.append(
                    Finding(
                        Severity.ERROR,
                        "NEGATIVE_VALUES_IN_LINEAR_SCALE",
                        "matrix",
                        None,
                        "values",
                        f"Negative value ({min_val:.4f}) detected in declared linear scale ({provenance.expression_unit.value}). Linear expression values must be non-negative.",
                    )
                )

        # Offset log scale non-negativity check (log2(x + c) for x >= 0, c >= 1)
        if provenance.transformation_scale == TransformationScale.LOG2 and provenance.pseudocount is not None and provenance.pseudocount >= 1.0:
            if min_val < 0.0:
                findings.append(
                    Finding(
                        Severity.ERROR,
                        "NEGATIVE_VALUES_IN_OFFSET_LOG_SCALE",
                        "matrix",
                        None,
                        "values",
                        f"Negative value ({min_val:.4f}) detected in declared log2(x+{provenance.pseudocount}) scale. Values must be non-negative when pseudocount >= 1.",
                    )
                )

        # Suspicious scale range warnings (non-blocking, preserving declared scale)
        if provenance.transformation_scale == TransformationScale.LINEAR and matrix.number_of_genes >= 1000:
            if max_val < 25.0 and provenance.expression_unit in {ExpressionUnit.TPM, ExpressionUnit.FPKM, ExpressionUnit.CPM}:
                findings.append(
                    Finding(
                        Severity.WARNING,
                        "SUSPICIOUS_EXPRESSION_RANGE",
                        "matrix",
                        None,
                        "values",
                        f"Maximum observed value ({max_val:.2f}) is unusually low for linear {provenance.expression_unit.value}, but declared linear scale will be respected without guessing.",
                    )
                )
        if provenance.is_log_transformed:
            if max_val > 1000.0:
                findings.append(
                    Finding(
                        Severity.WARNING,
                        "SUSPICIOUS_EXPRESSION_RANGE",
                        "matrix",
                        None,
                        "values",
                        f"Maximum observed value ({max_val:.2f}) is unusually high for log-transformed expression.",
                    )
                )

    ordered_findings_tuple = ordered(findings)
    overall_status = status_for(ordered_findings_tuple)
    is_eligible = overall_status != Status.FAIL

    return ProvenanceGateResult(
        status=overall_status,
        findings=ordered_findings_tuple,
        provenance=provenance,
        is_inferentially_eligible=is_eligible,
        deseq2_eligible=False,
        deseq2_ineligibility_reason="DESeq2 GLM requires raw unnormalized integer counts. Continuous processed expression matrices cannot be converted to integer counts or fed into DESeq2.",
        limma_trend_eligible=is_eligible,
        target_analysis_engine="limma_trend" if is_eligible else None,
        blocked_analysis_engines=("DESeq2",),
    )


def _load_reference(
    reference: GeneReference | None = None,
    config: dict[str, Any] | None = None,
) -> GeneReference:
    """Load canonical human gene reference catalog."""
    if reference is not None:
        return reference
    try:
        from src.gene_harmonization.config import load_config
        cfg = config or load_config()
        return GeneReference.load_from_config(cfg)
    except Exception:
        fixture_path = Path("tests/fixtures/ensembl_human_genes_test_fixture.tsv")
        if fixture_path.exists():
            cfg = {
                "reference": {
                    "file_path": str(fixture_path),
                    "reference_id": "ensembl_human_genes_test_fixture",
                    "reference_version": "112_fixture",
                    "genome_build": "GRCh38.p14",
                    "source_authority": "Ensembl / EMBL-EBI",
                }
            }
            return GeneReference.load_from_config(cfg)
        raise


def ingest_processed_expression_matrix(
    matrix_path: str | Path,
    metadata: dict[str, Any],
    sheet_name: str | None = None,
    asset_id: str | None = None,
    assay_id: str | None = None,
    provenance_overrides: dict[str, Any] | None = None,
    gene_reference: GeneReference | None = None,
    config: dict[str, Any] | None = None,
) -> ProcessedExpressionDataset:
    """Ingest, validate, QC, and provenance-gate a processed continuous expression matrix.

    Produces a canonical ProcessedExpressionDataset that establishes processing provenance,
    verifies scientific compatibility, strictly prevents accidental DESeq2 routing,
    and preserves original floating-point values for the future U3.3 analysis engine.
    """
    path_p = Path(matrix_path)

    # 1. Intake
    intake_data = read_processed_matrix(path_p, sheet_name=sheet_name)
    matrix = intake_data.matrix

    # 2. Quality Control
    qc_res = validate_processed_matrix(
        matrix_input=matrix if matrix is not None else path_p,
        metadata=metadata,
        asset_id=asset_id,
        assay_id=assay_id,
        config=config,
    )

    # 3. Provenance Gate
    gate_res = evaluate_provenance_gate(
        matrix=matrix,
        metadata=metadata,
        asset_id=asset_id,
        provenance_overrides=provenance_overrides,
        qc_result=qc_res,
    )

    # 4. Gene Harmonization
    harmonized_list: list[HarmonizedGene] = []
    if matrix is not None:
        try:
            ref = _load_reference(reference=gene_reference, config=config)
            for idx, g_id in enumerate(matrix.gene_ids, start=1):
                raw_id = g_id.strip()
                if raw_id.startswith("ENSG"):
                    id_type = GeneIdentifierType.ENSEMBL_GENE_ID
                elif raw_id.isdigit():
                    id_type = GeneIdentifierType.ENTREZ_GENE_ID
                elif raw_id.upper().startswith("HGNC:"):
                    id_type = GeneIdentifierType.HGNC_ID
                else:
                    id_type = GeneIdentifierType.GENE_SYMBOL

                lookup = ref.lookup(raw_id, id_type)

                harmonized_list.append(
                    HarmonizedGene(
                        source_index=idx,
                        original_gene_id=raw_id,
                        identifier_type=id_type.value,
                        canonical_gene_id=lookup.canonical_gene_id,
                        approved_symbol=lookup.approved_symbol,
                        mapping_status=lookup.status,
                        mapping_reason=lookup.reason,
                        candidate_canonical_ids=lookup.candidate_canonical_ids,
                        reference_id=ref.reference_id,
                        reference_version=ref.reference_version,
                        source_asset_id=asset_id or path_p.name,
                    )
                )
        except Exception:
            # Fallback if gene reference cannot be loaded: preserve original identifiers as unmapped
            for idx, g_id in enumerate(matrix.gene_ids, start=1):
                harmonized_list.append(
                    HarmonizedGene(
                        source_index=idx,
                        original_gene_id=g_id,
                        identifier_type="unknown",
                        canonical_gene_id=None,
                        approved_symbol=None,
                        mapping_status=GeneMappingStatus.UNMAPPED,
                        mapping_reason="reference_unavailable",
                        candidate_canonical_ids=(),
                        reference_id="none",
                        reference_version="none",
                        source_asset_id=asset_id or path_p.name,
                    )
                )

    harmonized_genes = tuple(harmonized_list)
    total_genes = matrix.number_of_genes if matrix else 0
    eligible_genes = sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED)
    excluded_genes = total_genes - eligible_genes

    # Combined status
    combined_findings = list(qc_res.findings) + list(gate_res.findings)
    overall_status = status_for(ordered(combined_findings))
    is_inferentially_eligible = qc_res.is_valid and gate_res.is_inferentially_eligible

    study_id = metadata.get("studies", [{}])[0].get("study_id") or metadata.get("study_id") or "processed_study"

    return ProcessedExpressionDataset(
        dataset_id=f"dataset_{study_id}",
        matrix=matrix if matrix is not None else ProcessedExpressionMatrix((), (), ()),
        provenance=gate_res.provenance,
        qc_result=qc_res,
        gate_result=gate_res,
        harmonized_genes=harmonized_genes,
        eligible_canonical_genes_count=eligible_genes,
        excluded_genes_count=excluded_genes,
        total_input_genes_count=total_genes,
        status=overall_status,
        is_inferentially_eligible=is_inferentially_eligible,
        target_analysis_engine=gate_res.target_analysis_engine if is_inferentially_eligible else None,
        blocked_analysis_engines=("DESeq2",),
    )
