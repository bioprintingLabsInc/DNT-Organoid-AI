"""Downstream pipeline connection bridge for PRJNA802072 validated counts.

Connects validated raw gene counts to:
- Step 6A: Quality Control (src.qc.validator.validate_bulk_counts)
- Step 6B: Gene Harmonization (src.gene_harmonization.harmonizer.harmonize_bulk_counts)
- Step 6C: Normalization (src.normalization.normalizer.normalize_bulk_dataset)
- Step 7B: Differential Expression (src.differential_expression.pipeline.DESeq2DifferentialExpressionPipeline)
- Stage 8C: Pathway & Signaling Analysis (src.features.pathway.pipeline.PathwayAnalysisPipeline)

Governed under verified biological replication criteria:
- Technical pilot (N=1 vs N=1) completes QC, Harmonization, and Normalization, but pauses before DESeq2
  because GLM dispersion estimation requires verified independent biological replicates (n >= 2).
- Full 12-library cohort (N=6 vs N=6) executes full end-to-end differential expression and pathway analysis.
"""

from copy import deepcopy
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import pandas as pd

from .manifest import APPROVED_PILOT_RUNS, RUN_TO_SAMPLE_MAP, TARGET_BIOPROJECT


@dataclass(frozen=True, slots=True)
class DownstreamBridgeResult:
    status: str
    qc_passed: bool
    qc_metrics: dict[str, Any]
    harmonization_metrics: dict[str, Any]
    normalization_metrics: dict[str, Any]
    deseq2_status: str
    deseq2_metrics: dict[str, Any]
    pathway_status: str
    pathway_metrics: dict[str, Any]
    biological_replicate_warning: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _reported(value: Any) -> dict[str, Any]:
    return {
        "original_value": value,
        "normalized_value": value,
        "value_status": "reported",
        "assertion_provenance_ids": [],
    }


def _attach_provenance(document: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    result = deepcopy(document)
    result["provenance"] = []
    identifiers = {
        "studies": "study_id",
        "experiments": "experiment_id",
        "conditions": "condition_id",
        "biological_sources": "biological_source_id",
        "organoid_contexts": "organoid_context_id",
        "samples": "sample_id",
        "exposures": "exposure_id",
        "treatment_control_relationships": "relationship_id",
        "sample_comparison_exceptions": "exception_id",
        "sequencing_assays": "assay_id",
        "input_data_assets": "asset_id",
        "dnt_reference_evidence": "evidence_id",
    }

    def attach(value: Any, entity: str, entity_id: str, field: str) -> None:
        if isinstance(value, dict):
            if "value_status" in value:
                prov_id = f"prov_{len(result['provenance']) + 1}"
                value["assertion_provenance_ids"] = [prov_id]
                result["provenance"].append({
                    "provenance_id": prov_id,
                    "entity_type": entity,
                    "entity_id": entity_id,
                    "field_path": field,
                    "source_kind": "source_file",
                    "source_reference": "canonical_metadata.yaml",
                })
            else:
                for nested in value.values():
                    attach(nested, entity, entity_id, field)
        elif isinstance(value, list):
            for nested in value:
                attach(nested, entity, entity_id, field)

    for entity, id_field in identifiers.items():
        if entity in result:
            for record in result[entity]:
                for field, value in record.items():
                    attach(value, entity, record.get(id_field, ""), field)
    return result


def build_canonical_metadata_for_runs(
    runs: tuple[str, ...],
    study_id: str = TARGET_BIOPROJECT,
    platform: str = "Illumina NovaSeq 6000",
) -> dict[str, Any]:
    """Construct canonical ingestion metadata structure for approved PRJNA802072 runs."""
    sample_records = []
    for r in runs:
        meta = RUN_TO_SAMPLE_MAP.get(r)
        is_ctrl = (meta.condition_group == "control") if meta else (r == "SRR17834920")
        cond_id = "control" if is_ctrl else "treatment_250um"
        dose = meta.vpa_dose_um if meta else (0.0 if is_ctrl else 250.0)
        day = meta.developmental_day if meta else 28

        sample_records.append({
            "sample_id": r,
            "experiment_id": f"{study_id}_exp1",
            "condition_id": cond_id,
            "biological_source_id": "line_1",
            "organoid_context_id": "cortical_1",
            "sample_alias": meta.sample_title if meta else r,
            "biosample_accession": meta.biosample_accession if meta else f"SAMN_{r}",
            "species": _reported("Homo sapiens"),
            "tissue_or_cell_type": _reported("Human Brain Cortical Organoid"),
            "developmental_day": _reported(day),
            "biological_replicate_id": _reported(f"rep_{r}"),
            "technical_replicate_id": None,
            "treatment_group": "control" if is_ctrl else "treatment",
            "compound_name": _reported("Valproic Acid"),
            "compound_concentration_um": _reported(dose),
            "exposure_duration_hours": _reported(72.0),
            "metadata_validation_status": "not_validated",
        })

    assay_id = f"{study_id}_assay"
    asset_id = f"{study_id}_counts"

    raw_doc = {
        "studies": [{
            "study_id": study_id,
            "source_type": "public",
            "metadata_validation_status": "not_validated",
        }],
        "experiments": [{
            "experiment_id": f"{study_id}_exp1",
            "study_id": study_id,
            "metadata_validation_status": "not_validated",
        }],
        "conditions": [
            {
                "condition_id": "control",
                "experiment_id": f"{study_id}_exp1",
                "condition_label": _reported("Vehicle Control (0 uM VPA)"),
                "treatment_control_status": _reported("control"),
                "control_type": _reported("vehicle"),
                "metadata_validation_status": "not_validated",
            },
            {
                "condition_id": "treatment_250um",
                "experiment_id": f"{study_id}_exp1",
                "condition_label": _reported("Valproic Acid 250 uM"),
                "treatment_control_status": _reported("treatment"),
                "metadata_validation_status": "not_validated",
            },
        ],
        "biological_sources": [{
            "biological_source_id": "line_1",
            "species": _reported("Homo sapiens"),
            "metadata_validation_status": "not_validated",
        }],
        "organoid_contexts": [{
            "organoid_context_id": "cortical_1",
            "organoid_type": _reported("cerebral organoid"),
            "brain_region_or_model_identity": _reported("cortical organoid"),
            "metadata_validation_status": "not_validated",
        }],
        "samples": sample_records,
        "exposures": [
            {
                "exposure_id": "control_exposure",
                "condition_id": "control",
                "exposure_scope": "condition_planned",
                "vehicle": _reported("Vehicle"),
                "developmental_age_or_stage_at_exposure": _reported("day 25"),
                "metadata_validation_status": "not_validated",
            },
            {
                "exposure_id": "treatment_exposure",
                "condition_id": "treatment_250um",
                "exposure_scope": "condition_planned",
                "agent_name": _reported("Valproic Acid"),
                "agent_identifier": _reported("pubchem:3121"),
                "vehicle": _reported("Vehicle"),
                "concentration_or_dose": _reported(250.0),
                "concentration_or_dose_unit": _reported("µM"),
                "exposure_start_time_or_stage": _reported("day 25"),
                "developmental_age_or_stage_at_exposure": _reported("day 25"),
                "exposure_duration": _reported(72),
                "exposure_duration_unit": _reported("hour"),
                "metadata_validation_status": "not_validated",
            },
        ],
        "treatment_control_relationships": [{
            "relationship_id": "comparison_1",
            "treatment_condition_ids": ["treatment_250um"],
            "matched_control_condition_ids": ["control"],
            "matching_basis": _reported("same organoid line, day 28 collection, 72h exposure"),
            "metadata_validation_status": "not_validated",
        }],
        "sample_comparison_exceptions": [],
        "sequencing_assays": [{
            "assay_id": assay_id,
            "sample_ids": list(runs),
            "rna_seq_modality": _reported("bulk_rna_seq"),
            "metadata_validation_status": "not_validated",
        }],
        "input_data_assets": [{
            "asset_id": asset_id,
            "assay_ids": [assay_id],
            "input_data_format": _reported("raw_counts"),
            "gene_identifier_type": _reported("Ensembl Gene ID"),
            "source_file_identifier_or_reference": _reported("counts.tsv"),
            "availability": "available",
            "metadata_validation_status": "not_validated",
        }],
        "provenance": [],
        "dnt_reference_evidence": [],
    }
    return _attach_provenance(raw_doc)


def execute_downstream_pipeline(
    counts_matrix_path: str | Path,
    runs: tuple[str, ...],
    study_id: str = TARGET_BIOPROJECT,
    contrast_name: str = "Day28_VPA_250uM_vs_Control",
    output_dir: str | Path | None = None,
) -> DownstreamBridgeResult:
    """Execute validated count matrix through existing frozen scientific modules.

    Applies strict biological replication guard:
    - If either treatment or control group has n < 2 samples, DESeq2 is blocked
      to prevent invalid dispersion estimates.
    """
    matrix_p = Path(counts_matrix_path).resolve()
    if not matrix_p.is_file():
        raise FileNotFoundError(f"Count matrix file not found: {matrix_p}")

    out_p = Path(output_dir).resolve() if output_dir else matrix_p.parent / "downstream_analysis"
    out_p.mkdir(parents=True, exist_ok=True)

    metadata = build_canonical_metadata_for_runs(runs, study_id=study_id)
    assay_id = f"{study_id}_assay"
    asset_id = f"{study_id}_counts"

    # -------------------------------------------------------------------------
    # Step 6A: Quality Control & Matrix Integrity Check
    # -------------------------------------------------------------------------
    from src.qc.validator import validate_bulk_counts
    from src.qc.result import Status as QCStatus

    qc_res = validate_bulk_counts(
        path=matrix_p,
        metadata=metadata,
        asset_id=asset_id,
        assay_id=assay_id,
    )
    qc_passed = qc_res.status != QCStatus.FAIL
    qc_metrics = {
        "status": str(qc_res.status),
        "total_genes": qc_res.number_of_genes,
        "total_samples": qc_res.number_of_samples,
        "findings_count": len(qc_res.findings),
    }

    # -------------------------------------------------------------------------
    # Step 6B: Gene Harmonization
    # -------------------------------------------------------------------------
    from src.gene_harmonization.harmonizer import harmonize_bulk_counts

    harm_metrics = {
        "status": "HARMONIZATION_INITIALIZED",
        "input_gene_count": 0,
        "mapped_genes_count": 0,
        "coverage_fraction": 0.0,
    }
    harm_res = None
    try:
        harm_res = harmonize_bulk_counts(
            matrix=matrix_p,
            metadata=metadata,
            asset_id=asset_id,
            assay_id=assay_id,
            qc_result=qc_res,
            allow_needs_review_qc=True,
        )
        harm_metrics["status"] = str(harm_res.status)
        if hasattr(harm_res, "summary"):
            harm_metrics["input_gene_count"] = harm_res.summary.total_input_genes
            harm_metrics["mapped_genes_count"] = harm_res.summary.uniquely_mapped_genes
            harm_metrics["coverage_fraction"] = harm_res.summary.canonical_coverage
    except Exception as e:
        harm_metrics["status"] = f"HARMONIZATION_DEFERRED: {e}"

    # -------------------------------------------------------------------------
    # Step 6C: Bulk Normalization
    # -------------------------------------------------------------------------
    from src.normalization.normalizer import normalize_bulk_dataset

    norm_metrics = {
        "status": "NORMALIZATION_INITIALIZED",
        "method": "median_of_ratios",
    }
    if harm_res is not None:
        try:
            norm_res = normalize_bulk_dataset(
                harmonized=harm_res,
                metadata=metadata,
            )
            norm_metrics["status"] = str(norm_res.status)
        except Exception as e:
            norm_metrics["status"] = f"NORMALIZATION_DEFERRED: {e}"

    # -------------------------------------------------------------------------
    # Step 7B: Differential Expression (DESeq2)
    # -------------------------------------------------------------------------
    controls = [r for r in runs if RUN_TO_SAMPLE_MAP.get(r, None) and RUN_TO_SAMPLE_MAP[r].condition_group == "control"]
    treated = [r for r in runs if RUN_TO_SAMPLE_MAP.get(r, None) and RUN_TO_SAMPLE_MAP[r].condition_group == "treatment"]

    n_ctrl = len(controls)
    n_trt = len(treated)

    deseq2_status = "NOT_EVALUATED"
    deseq2_metrics: dict[str, Any] = {}
    pathway_status = "NOT_EVALUATED"
    pathway_metrics: dict[str, Any] = {}
    replicate_warning = None

    if n_ctrl < 2 or n_trt < 2:
        deseq2_status = "BLOCKED_INSUFFICIENT_BIOLOGICAL_REPLICATION"
        replicate_warning = (
            f"DESeq2 GLM dispersion estimation requires n >= 2 biological replicates per group. "
            f"Observed: controls={n_ctrl}, treated={n_trt}. Execution safely deferred until full cohort."
        )
        deseq2_metrics = {
            "status": deseq2_status,
            "controls_count": n_ctrl,
            "treated_count": n_trt,
            "reason": replicate_warning,
        }
    else:
        # Full cohort (n >= 2)
        from src.differential_expression.pipeline import DESeq2DifferentialExpressionPipeline
        from src.differential_expression.design import ContrastSpecification

        diff_pipeline = DESeq2DifferentialExpressionPipeline()
        spec = ContrastSpecification(
            contrast_id=contrast_name,
            test_condition="treatment",
            baseline_condition="control",
            design_formula="~ condition",
        )
        try:
            diff_res = diff_pipeline.run_contrast(
                counts_path=matrix_p,
                metadata=metadata,
                contrast_spec=spec,
            )
            deseq2_status = "DESEQ2_COMPLETED"
            deseq2_metrics = {
                "status": "DESEQ2_SUCCESS",
                "contrast_id": contrast_name,
                "tested_genes": getattr(diff_res, "total_tested_genes", 0),
                "significant_degs": getattr(diff_res, "significant_genes_count", 0),
            }

            # -----------------------------------------------------------------
            # Stage 8C: Pathway & Signaling Analysis
            # -----------------------------------------------------------------
            from src.features.pathway.pipeline import PathwayAnalysisPipeline

            path_pipeline = PathwayAnalysisPipeline()
            path_res = path_pipeline.run_gsea(diff_res)
            pathway_status = "PATHWAY_ANALYSIS_SUCCESS"
            pathway_metrics = {
                "status": "GSEA_SUCCESS",
                "go_bp_significant": getattr(path_res, "go_bp_significant_count", 0),
                "kegg_significant": getattr(path_res, "kegg_significant_count", 0),
            }
        except Exception as e:
            deseq2_status = f"DESEQ2_ERROR: {e}"
            deseq2_metrics = {"error": str(e)}

    return DownstreamBridgeResult(
        status="DOWNSTREAM_EVALUATION_COMPLETED",
        qc_passed=qc_passed,
        qc_metrics=qc_metrics,
        harmonization_metrics=harm_metrics,
        normalization_metrics=norm_metrics,
        deseq2_status=deseq2_status,
        deseq2_metrics=deseq2_metrics,
        pathway_status=pathway_status,
        pathway_metrics=pathway_metrics,
        biological_replicate_warning=replicate_warning,
    )
