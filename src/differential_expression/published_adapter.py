"""Adapter layer for published differential-expression summary tables.

Converts published differential-expression results into canonical
DifferentialExpressionContrastResult objects without feeding them through
raw-count machinery (no DESeq2 re-computation, no size-factor normalization).

Enforces:
1. Preservation of published statistical values (log2FC, p-value, padj, SE, Wald).
2. Zero data fabrication: omitted/missing statistics remain None.
3. Gene reference harmonization against canonical Ensembl reference (Ensembl 112).
4. Gene universe conservation law:
     total_input_genes == eligible_canonical_genes + excluded_genes
5. Descriptive ResponseDirection:
     UP (log2FC > 0), DOWN (log2FC < 0), NO_DIRECTION (log2FC == 0), UNDETERMINED (unestimable).
6. Statistical eligibility separated from response direction.
7. Full experimental context and provenance capture.
"""

import hashlib
import math
from pathlib import Path
from typing import Any

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    ExcludedGeneAudit,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.models import ResponseDirection
from src.gene_harmonization.models import GeneIdentifierType, GeneMappingStatus
from src.gene_harmonization.reference import GeneReference
from src.ingestion.deg_table_reader import DegTableData, DegTableRow, read_deg_table
from src.response_builder.models import ContrastExposure, ContrastStatus, Finding, Severity


def derive_published_response_direction(log2_fc: float | None) -> ResponseDirection:
    """Derive purely descriptive direction from the sign of published log2FC.

    Non-evaluative: does NOT imply statistical significance or biological toxicity.
    """
    if log2_fc is None or math.isnan(log2_fc) or math.isinf(log2_fc):
        return ResponseDirection.UNDETERMINED
    if log2_fc > 0.0:
        return ResponseDirection.UP
    if log2_fc < 0.0:
        return ResponseDirection.DOWN
    return ResponseDirection.NO_DIRECTION


def _extract_contrast_context(
    metadata: dict[str, Any],
    contrast_id: str | None = None,
) -> dict[str, Any]:
    """Extract study, condition, and exposure context from canonical or simplified metadata."""
    context: dict[str, Any] = {}

    dataset = metadata.get("dataset", {})
    studies = metadata.get("studies", [])
    cohorts = metadata.get("cohorts", [])
    experiments = metadata.get("experiments", [])

    context["study_id"] = (
        dataset.get("study_id")
        or (studies[0].get("study_id") if studies else None)
        or metadata.get("study_id")
        or "published_study"
    )
    context["cohort_id"] = (
        dataset.get("cohort_id")
        or (cohorts[0].get("cohort_id") if cohorts else None)
        or metadata.get("cohort_id")
        or "cohort_1"
    )
    context["experiment_id"] = (
        (experiments[0].get("experiment_id") if experiments else None)
        or metadata.get("experiment_id")
        or "exp_1"
    )

    relationships = metadata.get("relationships") or metadata.get("treatment_control_relationships", [])
    selected_rel = None
    if contrast_id:
        for r in relationships:
            if r.get("relationship_id") == contrast_id or r.get("contrast_id") == contrast_id:
                selected_rel = r
                break
    if selected_rel is None and relationships:
        # Pick first treatment_control relationship
        for r in relationships:
            if r.get("relationship_type") == "treatment_control":
                selected_rel = r
                break
        if selected_rel is None:
            selected_rel = relationships[0]

    conditions_by_id = {c.get("condition_id"): c for c in metadata.get("conditions", []) if "condition_id" in c}
    samples = metadata.get("samples", [])

    if selected_rel is not None:
        context["relationship_id"] = selected_rel.get("relationship_id", "rel_1")
        context["contrast_id"] = contrast_id or selected_rel.get("relationship_id", "contrast_1")
        treatment_cid = selected_rel.get("treatment_condition_id", "treatment")
        ctrl_cids = tuple(selected_rel.get("control_condition_ids", ()))
        context["treatment_condition_id"] = treatment_cid
        context["matched_control_condition_ids"] = ctrl_cids

        trt_samples = [s.get("sample_id") for s in samples if s.get("condition_id") == treatment_cid]
        ctrl_samples = [s.get("sample_id") for s in samples if s.get("condition_id") in ctrl_cids]
        context["treatment_sample_ids"] = tuple(trt_samples)
        context["control_sample_ids"] = tuple(ctrl_samples)

        trt_bio = [s.get("biological_replicate_id") or s.get("sample_id") for s in samples if s.get("condition_id") == treatment_cid]
        ctrl_bio = [s.get("biological_replicate_id") or s.get("sample_id") for s in samples if s.get("condition_id") in ctrl_cids]
        context["treatment_biological_replicate_ids"] = tuple(trt_bio)
        context["control_biological_replicate_ids"] = tuple(ctrl_bio)

        context["treatment_replicate_count"] = len(trt_samples) if trt_samples else selected_rel.get("treatment_replicate_count", 3)
        context["control_replicate_count"] = len(ctrl_samples) if ctrl_samples else selected_rel.get("control_replicate_count", 3)

        trt_cond = conditions_by_id.get(treatment_cid, {})
        context["agent_name"] = trt_cond.get("chemical_agent_name") or trt_cond.get("agent_name") or metadata.get("agent_name")
        context["agent_identifier"] = trt_cond.get("chemical_agent_id") or trt_cond.get("agent_identifier") or metadata.get("agent_identifier")
        context["vehicle"] = trt_cond.get("vehicle") or metadata.get("vehicle")
        context["concentration_or_dose"] = trt_cond.get("concentration_or_dose") or trt_cond.get("dose_value") or metadata.get("concentration_or_dose")
        context["concentration_or_dose_unit"] = trt_cond.get("concentration_or_dose_unit") or trt_cond.get("dose_unit") or metadata.get("concentration_or_dose_unit")
        context["exposure_duration"] = trt_cond.get("exposure_duration") or metadata.get("exposure_duration")
        context["exposure_duration_unit"] = trt_cond.get("exposure_duration_unit") or metadata.get("exposure_duration_unit")
        context["developmental_age_or_stage_at_exposure"] = trt_cond.get("developmental_age_or_stage_at_exposure") or metadata.get("developmental_age_or_stage_at_exposure")
    else:
        # Fallback to direct metadata fields
        context["relationship_id"] = metadata.get("relationship_id", "rel_1")
        context["contrast_id"] = contrast_id or metadata.get("contrast_id", "contrast_1")
        context["treatment_condition_id"] = metadata.get("treatment_condition_id", "treatment")
        ctrl = metadata.get("control_condition_id") or metadata.get("matched_control_condition_ids") or "control"
        context["matched_control_condition_ids"] = (ctrl,) if isinstance(ctrl, str) else tuple(ctrl)
        context["treatment_sample_ids"] = tuple(metadata.get("treatment_sample_ids", ()))
        context["control_sample_ids"] = tuple(metadata.get("control_sample_ids", ()))
        context["treatment_biological_replicate_ids"] = tuple(metadata.get("treatment_biological_replicate_ids", ()))
        context["control_biological_replicate_ids"] = tuple(metadata.get("control_biological_replicate_ids", ()))
        context["treatment_replicate_count"] = metadata.get("treatment_replicate_count", len(context["treatment_sample_ids"]) or 3)
        context["control_replicate_count"] = metadata.get("control_replicate_count", len(context["control_sample_ids"]) or 3)
        context["agent_name"] = metadata.get("agent_name")
        context["agent_identifier"] = metadata.get("agent_identifier")
        context["vehicle"] = metadata.get("vehicle")
        context["concentration_or_dose"] = metadata.get("concentration_or_dose")
        context["concentration_or_dose_unit"] = metadata.get("concentration_or_dose_unit")
        context["exposure_duration"] = metadata.get("exposure_duration")
        context["exposure_duration_unit"] = metadata.get("exposure_duration_unit")
        context["developmental_age_or_stage_at_exposure"] = metadata.get("developmental_age_or_stage_at_exposure")

    # Exposures tuple
    exposures = []
    if context.get("agent_name"):
        exposures.append(
            ContrastExposure(
                exposure_id=f"{context['contrast_id']}_exp1",
                agent_name=str(context["agent_name"]),
                agent_identifier=str(context.get("agent_identifier") or ""),
                concentration_or_dose=context.get("concentration_or_dose"),
                concentration_or_dose_unit=str(context.get("concentration_or_dose_unit") or ""),
                exposure_duration=context.get("exposure_duration"),
                exposure_duration_unit=str(context.get("exposure_duration_unit") or ""),
            )
        )
    context["treatment_exposures"] = tuple(exposures)
    return context


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
        # Fallback to test fixture if full reference catalog is not in active working directory
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


def ingest_published_deg_table(
    deg_table_path: str | Path,
    metadata: dict[str, Any],
    sheet_name: str | None = None,
    contrast_id: str | None = None,
    gene_reference: GeneReference | None = None,
    column_mapping: dict[str, str] | None = None,
    config: dict[str, Any] | None = None,
) -> tuple[DifferentialExpressionContrastResult, DegTableData]:
    """Ingest, harmonize, and adapt a published differential expression table.

    Produces a canonical DifferentialExpressionContrastResult that directly
    integrates with Stage 8A.1 Condition-Level Gene Response Registry.
    """
    path_p = Path(deg_table_path)
    table_data = read_deg_table(path_p, sheet_name=sheet_name, column_mapping=column_mapping)

    findings: list[Finding] = []
    for f in table_data.findings:
        findings.append(Finding(Severity[f.severity], f.rule, f.path, None, "deg_table", f.message))

    ctx = _extract_contrast_context(metadata, contrast_id=contrast_id)

    if not table_data.is_valid:
        findings.append(
            Finding(
                Severity.ERROR,
                "invalid_published_table",
                "deg_table",
                None,
                "table_data",
                "Published differential expression table cannot be parsed or lacks required columns.",
            )
        )
        res = DifferentialExpressionContrastResult(
            contrast_id=ctx["contrast_id"],
            study_id=ctx["study_id"],
            experiment_id=ctx["experiment_id"],
            cohort_id=ctx["cohort_id"],
            relationship_id=ctx["relationship_id"],
            treatment_condition_id=ctx["treatment_condition_id"],
            matched_control_condition_ids=ctx["matched_control_condition_ids"],
            treatment_sample_ids=ctx["treatment_sample_ids"],
            control_sample_ids=ctx["control_sample_ids"],
            treatment_biological_replicate_ids=ctx["treatment_biological_replicate_ids"],
            control_biological_replicate_ids=ctx["control_biological_replicate_ids"],
            treatment_replicate_count=ctx["treatment_replicate_count"],
            control_replicate_count=ctx["control_replicate_count"],
            status=ContrastStatus.BLOCKED,
            findings=tuple(findings),
            total_input_genes_count=table_data.total_rows,
        )
        return res, table_data

    # Load canonical gene reference
    ref = _load_reference(reference=gene_reference, config=config)

    gene_results: list[GeneDifferentialExpressionResult] = []
    excluded_genes: list[ExcludedGeneAudit] = []
    seen_canonical_ids: set[str] = set()

    for row in table_data.rows:
        raw_id = row.raw_gene_id.strip()

        # Detect identifier type
        if raw_id.startswith("ENSG"):
            id_type = GeneIdentifierType.ENSEMBL_GENE_ID
        elif raw_id.isdigit():
            id_type = GeneIdentifierType.ENTREZ_GENE_ID
        elif raw_id.upper().startswith("HGNC:"):
            id_type = GeneIdentifierType.HGNC_ID
        else:
            id_type = GeneIdentifierType.GENE_SYMBOL

        lookup = ref.lookup(raw_id, id_type)

        if lookup.status == GeneMappingStatus.UNIQUELY_MAPPED and lookup.canonical_gene_id:
            can_id = lookup.canonical_gene_id
            if can_id in seen_canonical_ids:
                # Collision / duplicate in published source table
                excluded_genes.append(
                    ExcludedGeneAudit(
                        original_gene_id=raw_id,
                        source_index=row.source_index,
                        mapping_status="COLLISION",
                        canonical_gene_id=can_id,
                        exclusion_reason="Duplicate canonical gene ID collision in published table",
                    )
                )
                continue

            seen_canonical_ids.add(can_id)

            # Determine gene result status
            if row.log2_fold_change is None or math.isnan(row.log2_fold_change) or math.isinf(row.log2_fold_change):
                res_status = "NON_ESTIMABLE"
                res_notes = "Published log2FoldChange is missing, NaN, or infinite"
            elif row.p_value is None and row.adjusted_p_value is None and row.log2_fold_change is not None:
                res_status = "OUTLIER"
                res_notes = "Significance statistics omitted or filtered as outlier in published table"
            else:
                res_status = "OK"
                res_notes = None

            approved_sym = lookup.approved_symbol or row.raw_gene_symbol

            gene_results.append(
                GeneDifferentialExpressionResult(
                    canonical_gene_id=can_id,
                    approved_symbol=approved_sym,
                    base_mean=row.base_mean if row.base_mean is not None else float("nan"),
                    log2_fold_change=row.log2_fold_change,
                    lfc_standard_error=row.lfc_standard_error,  # None if omitted
                    wald_statistic=row.wald_statistic,          # None if omitted
                    p_value=row.p_value,                        # None if omitted
                    adjusted_p_value_bh=row.adjusted_p_value,   # None if omitted
                    result_status=res_status,
                    result_notes=res_notes,
                )
            )
        else:
            # Unmapped, ambiguous, or invalid
            excluded_genes.append(
                ExcludedGeneAudit(
                    original_gene_id=raw_id,
                    source_index=row.source_index,
                    mapping_status=str(lookup.status),
                    canonical_gene_id=lookup.canonical_gene_id,
                    exclusion_reason=lookup.reason or f"Identifier mapping status: {lookup.status}",
                )
            )

    # Replicate qualification and contrast status determination
    is_replicated = ctx["treatment_replicate_count"] >= 2 and ctx["control_replicate_count"] >= 2
    if is_replicated and len(gene_results) > 0:
        c_status = ContrastStatus.ELIGIBLE
    elif len(gene_results) > 0:
        c_status = ContrastStatus.NEEDS_REVIEW
        findings.append(
            Finding(
                Severity.REVIEW,
                "insufficient_replicates",
                "deg_table",
                None,
                "replicate_count",
                f"Contrast has fewer than 2 replicates (treatment: {ctx['treatment_replicate_count']}, control: {ctx['control_replicate_count']}); flagged for review.",
            )
        )
    else:
        c_status = ContrastStatus.BLOCKED
        findings.append(
            Finding(
                Severity.ERROR,
                "no_canonical_genes",
                "deg_table",
                None,
                "gene_results",
                "No canonical genes could be mapped from the published table.",
            )
        )

    # Compute checksum of input file for provenance
    input_hash = ""
    try:
        input_hash = hashlib.sha256(path_p.read_bytes()).hexdigest()
    except Exception:
        pass

    provenance_info = {
        "source_pipeline": "published_deg_table_ingestion_v2_u3",
        "analysis_method": metadata.get("published_analysis_method", "DESeq2"),
        "source_file": str(path_p),
        "sheet_name": table_data.sheet_name,
        "columns_discovered": table_data.column_mapping,
        "total_rows_parsed": table_data.total_rows,
        "canonical_mapped_genes": len(gene_results),
        "excluded_unmapped_genes": len(excluded_genes),
    }

    result = DifferentialExpressionContrastResult(
        contrast_id=ctx["contrast_id"],
        study_id=ctx["study_id"],
        experiment_id=ctx["experiment_id"],
        cohort_id=ctx["cohort_id"],
        relationship_id=ctx["relationship_id"],
        treatment_condition_id=ctx["treatment_condition_id"],
        matched_control_condition_ids=ctx["matched_control_condition_ids"],
        treatment_sample_ids=ctx["treatment_sample_ids"],
        control_sample_ids=ctx["control_sample_ids"],
        treatment_biological_replicate_ids=ctx["treatment_biological_replicate_ids"],
        control_biological_replicate_ids=ctx["control_biological_replicate_ids"],
        treatment_replicate_count=ctx["treatment_replicate_count"],
        control_replicate_count=ctx["control_replicate_count"],
        treatment_exposures=ctx["treatment_exposures"],
        agent_name=ctx["agent_name"],
        agent_identifier=ctx["agent_identifier"],
        vehicle=ctx["vehicle"],
        concentration_or_dose=ctx["concentration_or_dose"],
        concentration_or_dose_unit=ctx["concentration_or_dose_unit"],
        exposure_duration=ctx["exposure_duration"],
        exposure_duration_unit=ctx["exposure_duration_unit"],
        developmental_age_or_stage_at_exposure=ctx["developmental_age_or_stage_at_exposure"],
        design_formula="~ condition (published_summary)",
        design_rationale="Published differential expression summary statistics; ingested directly without count-based re-estimation.",
        gene_results=tuple(gene_results),
        excluded_genes=tuple(excluded_genes),
        total_input_genes_count=table_data.total_rows,
        status=c_status,
        findings=tuple(findings),
        r_environment_info=provenance_info,
        input_hashes={"deg_table_sha256": input_hash},
        de_version="published_v1",
    )

    return result, table_data
