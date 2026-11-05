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
- Full 12-library cohort (N=6 vs N=6) completes QC, Harmonization, and Normalization, but inferential DESeq2
  and pathway analysis are gated (DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED) when differentiation
  batches and organoid pooling details are unverified, preventing ungrounded biological claims.
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
    metadata_validation_status: str | None = None
    metadata_validation_metrics: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _reported(value: Any) -> dict[str, Any]:
    return {
        "original_value": value,
        "normalized_value": value,
        "value_status": "reported",
        "assertion_provenance_ids": [],
    }


def _missing(status: str = "missing_not_reported") -> dict[str, Any]:
    return {
        "original_value": None,
        "normalized_value": None,
        "value_status": status,
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
    strict_schema: bool = False,
    verified_biological_replicates: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Construct canonical ingestion metadata structure for approved PRJNA802072 runs."""
    sample_records = []
    for r in runs:
        meta = RUN_TO_SAMPLE_MAP.get(r)
        is_ctrl = (meta.condition_group == "control") if meta else (r == "SRR17834920")
        cond_id = "control" if is_ctrl else "treatment_250um"
        dose = meta.vpa_dose_um if meta else (0.0 if is_ctrl else 250.0)
        day = meta.developmental_day if meta else 28

        # Biological replicate resolution
        if verified_biological_replicates and r in verified_biological_replicates:
            bio_rep = _reported(verified_biological_replicates[r])
        else:
            bio_rep = _missing("unknown")

        s_rec: dict[str, Any] = {
            "sample_id": r,
            "experiment_id": f"{study_id}_exp1",
            "condition_id": cond_id,
            "biological_source_id": "line_1",
            "organoid_context_id": "cortical_1",
            "biological_replicate_id": bio_rep,
            "technical_replicate_id": _missing("not_applicable"),
            "developmental_age_or_stage_at_collection": _reported("Day 28"),
            "collection_time_relative_to_exposure": _reported("288.0 hour"),
            "metadata_validation_status": "not_validated",
        }

        # Include legacy/convenience non-schema fields when strict_schema is False
        if not strict_schema:
            s_rec.update({
                "sample_alias": meta.sample_title if meta else r,
                "biosample_accession": meta.biosample_accession if meta else f"SAMN_{r}",
                "species": _reported("Homo sapiens"),
                "tissue_or_cell_type": _reported("Human Brain Cortical Organoid"),
                "developmental_day": _reported(day),
                "treatment_group": "control" if is_ctrl else "treatment",
                "compound_name": _reported("Vehicle") if is_ctrl else _reported("Valproic Acid"),
                "compound_concentration_um": _reported(dose),
                "exposure_duration_hours": _reported(288.0),
            })

        sample_records.append(s_rec)

    assay_id = f"{study_id}_assay"
    asset_id = f"{study_id}_counts"

    assay_dict: dict[str, Any] = {
        "assay_id": assay_id,
        "sample_ids": list(runs),
        "rna_seq_modality": _reported("bulk_rna_seq"),
        "read_layout": _reported("paired_end"),
        "strandedness": _reported("unstranded"),
        "sequencing_platform": _reported(platform),
        "metadata_validation_status": "not_validated",
    }
    if not strict_schema:
        assay_dict["library_layout"] = _reported("paired_end")

    raw_doc = {
        "studies": [{
            "study_id": study_id,
            "source_type": "public",
            "public_accessions": [{
                "namespace": "BioProject",
                "identifier": study_id,
                "url": f"https://www.ncbi.nlm.nih.gov/bioproject/{study_id}",
            }],
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
            "source_material_type": _reported("hESC"),
            "cell_line_or_source_identifier": _reported("H9"),
            "metadata_validation_status": "not_validated",
        }],
        "organoid_contexts": [{
            "organoid_context_id": "cortical_1",
            "organoid_type": _reported("cerebral organoid"),
            "brain_region_or_model_identity": _reported("dorsal forebrain organoids"),
            "metadata_validation_status": "not_validated",
        }],
        "samples": sample_records,
        "exposures": [
            {
                "exposure_id": "control_exposure",
                "condition_id": "control",
                "exposure_scope": "condition_planned",
                "vehicle": _reported("distilled water"),
                "exposure_start_time_or_stage": _reported("day 16"),
                "developmental_age_or_stage_at_exposure": _reported("day 16"),
                "exposure_duration": _reported(288.0),
                "exposure_duration_unit": _reported("hour"),
                "metadata_validation_status": "not_validated",
            },
            {
                "exposure_id": "treatment_exposure",
                "condition_id": "treatment_250um",
                "exposure_scope": "condition_planned",
                "agent_name": _reported("Valproic Acid"),
                "agent_identifier": _reported("pubchem:3121"),
                "vehicle": _reported("distilled water"),
                "concentration_or_dose": _reported(250.0),
                "concentration_or_dose_unit": _reported("µM"),
                "exposure_start_time_or_stage": _reported("day 16"),
                "developmental_age_or_stage_at_exposure": _reported("day 16"),
                "exposure_duration": _reported(288.0),
                "exposure_duration_unit": _reported("hour"),
                "metadata_validation_status": "not_validated",
            },
        ],
        "treatment_control_relationships": [{
            "relationship_id": "comparison_1",
            "treatment_condition_ids": ["treatment_250um"],
            "matched_control_condition_ids": ["control"],
            "matching_basis": _reported("H9 hESC dorsal forebrain organoids, Day 16-28 (12 days / 288h) VPA exposure vs vehicle"),
            "metadata_validation_status": "not_validated",
        }],
        "sample_comparison_exceptions": [],
        "sequencing_assays": [assay_dict],
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
    verified_biological_replicates: dict[str, str] | None = None,
    canonical_metadata: dict[str, Any] | None = None,
) -> DownstreamBridgeResult:
    """Execute validated count matrix through existing frozen scientific modules.

    Applies strict biological replication guard:
    - If either treatment or control group has n < 2 samples, DESeq2 is blocked
      to prevent invalid dispersion estimates.
    - If biological replicates are unverified (e.g., differentiation batches and
      organoid pooling details are unknown), DESeq2 inferential GLM and pathway
      analysis are gated (DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED) while
      QC, Harmonization, and Normalization pass.
    """
    matrix_p = Path(counts_matrix_path).resolve()
    if not matrix_p.is_file():
        raise FileNotFoundError(f"Count matrix file not found: {matrix_p}")

    out_p = Path(output_dir).resolve() if output_dir else matrix_p.parent / "downstream_analysis"
    out_p.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # Pre-Flight: Universal Metadata Validation and Experimental Integrity Gate
    # -------------------------------------------------------------------------
    from src.metadata.validator import MetadataValidator
    from src.metadata.integrity import MetadataIntegrityChecker
    from src.metadata.result import ValidationStatus

    # Build schema-conformant canonical metadata document for validator
    if canonical_metadata is not None:
        schema_metadata = canonical_metadata
        metadata = canonical_metadata
    else:
        schema_metadata = build_canonical_metadata_for_runs(
            runs,
            study_id=study_id,
            strict_schema=True,
            verified_biological_replicates=verified_biological_replicates,
        )
        metadata = build_canonical_metadata_for_runs(
            runs,
            study_id=study_id,
            strict_schema=False,
            verified_biological_replicates=verified_biological_replicates,
        )

    validator = MetadataValidator()
    val_res = validator.validate(schema_metadata)

    checker = MetadataIntegrityChecker(min_biological_replicates=2)
    integrity_res = checker.check(schema_metadata)

    val_is_valid = (val_res.status != ValidationStatus.INVALID)
    metadata_val_status = str(val_res.status.value)
    metadata_val_metrics = {
        "validator_status": str(val_res.status.value),
        "validator_is_valid": val_is_valid,
        "validator_findings_count": len(val_res.findings),
        "integrity_status": str(integrity_res.status.value),
        "integrity_is_valid": integrity_res.is_valid,
        "integrity_findings_count": len(integrity_res.findings),
        "replicate_summary": integrity_res.replicate_summary,
        "model_consistency_summary": integrity_res.model_consistency_summary,
    }

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

    advisory_flags: list[str] = []
    if "SRR17834908" in runs:
        advisory_flags.append(
            "SRR17834908: spearman_correlation=0.8860 below 0.90 threshold, PC1/PC2 separation observed; "
            "retained for exploration without automated exclusion"
        )

    qc_metrics = {
        "status": str(qc_res.status),
        "total_genes": qc_res.number_of_genes,
        "total_samples": qc_res.number_of_samples,
        "findings_count": len(qc_res.findings),
        "advisory_flags": advisory_flags,
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
        if hasattr(harm_res, "summary") and harm_res.summary is not None:
            harm_metrics["input_gene_count"] = getattr(harm_res.summary, "total_genes", 0)
            harm_metrics["mapped_genes_count"] = getattr(harm_res.summary, "uniquely_mapped_genes", 0)
            harm_metrics["coverage_fraction"] = getattr(harm_res.summary, "percentage_uniquely_mapped", 0.0)
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
    # Step 7B: Differential Expression (DESeq2) Biological Replication Gate
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
    final_status = "DOWNSTREAM_EVALUATION_COMPLETED"

    if not val_is_valid:
        deseq2_status = "BLOCKED_INVALID_CANONICAL_METADATA"
        pathway_status = "BLOCKED_INVALID_CANONICAL_METADATA"
        err_findings = [f.message for f in val_res.findings if f.severity.value in ("ERROR", "CRITICAL")]
        if not err_findings:
            err_findings = [f.message for f in val_res.findings]
        replicate_warning = (
            f"Canonical metadata validation failed ({val_res.status.value}): "
            + "; ".join(err_findings[:3])
        )
        deseq2_metrics = {
            "status": deseq2_status,
            "reason": replicate_warning,
            "findings": [f.to_dict() for f in val_res.findings],
        }
        pathway_metrics = {
            "status": pathway_status,
            "reason": replicate_warning,
        }
        final_status = "DOWNSTREAM_BLOCKED_INVALID_METADATA"
    elif n_ctrl < 2 or n_trt < 2:
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
        pathway_status = "NOT_EVALUATED"
        pathway_metrics = {"reason": replicate_warning}
        final_status = "DOWNSTREAM_PARTIAL_INSUFFICIENT_REPLICATION"
    else:
        # Check whether biological replication has been verified
        has_unverified_replicates = (
            not integrity_res.is_valid
            or any(
                s.get("biological_replicate_id", {}).get("value_status") != "reported"
                for s in schema_metadata.get("samples", [])
            )
            or (
                verified_biological_replicates is None
                and any(
                    getattr(RUN_TO_SAMPLE_MAP.get(r), "biological_replicate_unit", "") == "UNVERIFIED"
                    for r in runs
                )
            )
        )

        if has_unverified_replicates:
            deseq2_status = "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
            pathway_status = "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
            replicate_warning = (
                "Biological replication is UNVERIFIED for PRJNA802072: differentiation batch IDs "
                "and organoid pooling details are not reported. SRA BioSample uniqueness alone does "
                "not establish independent biological replication. DESeq2 inferential GLM and pathway "
                "analysis are safely gated to prevent ungrounded biological claims."
            )
            deseq2_metrics = {
                "status": deseq2_status,
                "controls_count": n_ctrl,
                "treated_count": n_trt,
                "reason": replicate_warning,
            }
            pathway_metrics = {
                "status": pathway_status,
                "reason": replicate_warning,
            }
            final_status = "DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED"
        else:
            # Full cohort with verified biological replicates (n >= 2)
            # -------------------------------------------------------------
            # Step 7A: Treatment-versus-Matched-Control Contrast Building
            # -------------------------------------------------------------
            from src.response_builder.builder import build_response_contrasts
            from src.response_builder.models import ContrastStatus

            if norm_res is not None and harm_res is not None:
                contrast_ds = build_response_contrasts(norm_res, schema_metadata)
                primary_contrast = next((c for c in contrast_ds.contrasts if c.contrast_id == contrast_name), None)
                if primary_contrast is None and contrast_ds.contrasts:
                    primary_contrast = contrast_ds.contrasts[0]

                if primary_contrast is None:
                    deseq2_status = "BLOCKED_NO_CONTRASTS_BUILT"
                    deseq2_metrics = {
                        "status": deseq2_status,
                        "reasons": [f.message for f in contrast_ds.findings],
                    }
                    pathway_status = "NOT_EVALUATED"
                    final_status = "DOWNSTREAM_BLOCKED_INFERENTIALLY"
                elif primary_contrast.status != ContrastStatus.ELIGIBLE or not primary_contrast.is_inferentially_eligible:
                    deseq2_status = f"BLOCKED_{primary_contrast.status.value}"
                    deseq2_metrics = {
                        "status": deseq2_status,
                        "contrast_id": primary_contrast.contrast_id,
                        "reasons": [f.message for f in primary_contrast.findings],
                    }
                    pathway_status = "NOT_EVALUATED"
                    final_status = "DOWNSTREAM_BLOCKED_INFERENTIALLY"
                else:
                    # ---------------------------------------------------------
                    # Step 7B: Differential Expression (DESeq2 Wald Test)
                    # ---------------------------------------------------------
                    from src.differential_expression.pipeline import run_contrast_differential_expression

                    try:
                        de_res = run_contrast_differential_expression(
                            contrast=primary_contrast,
                            normalized_dataset=norm_res,
                            harmonized_dataset=harm_res,
                            sample_metadata=schema_metadata.get("samples", []),
                        )
                        if de_res.status == ContrastStatus.ELIGIBLE:
                            deseq2_status = "DESEQ2_COMPLETED"
                            sig_count = sum(1 for g in de_res.gene_results if g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05)
                            deseq2_metrics = {
                                "status": "DESEQ2_SUCCESS",
                                "contrast_id": primary_contrast.contrast_id,
                                "tested_genes": len(de_res.gene_results),
                                "significant_degs": sig_count,
                                "design_formula": de_res.design_formula,
                            }

                            # -----------------------------------------------------
                            # Stage 8C: Pathway & Signaling Analysis
                            # -----------------------------------------------------
                            from src.features.pathway.pipeline import run_pathway_enrichment_for_contrast

                            try:
                                path_res = run_pathway_enrichment_for_contrast(de_res)
                                ranked_count = len(path_res.ranked_genes.ranked_genes)
                                go_count = len(path_res.go_results)
                                kegg_count = len(path_res.kegg_results)
                                sig_pathways = (
                                    path_res.qc.go_significant_pathways_count
                                    + path_res.qc.kegg_significant_pathways_count
                                )

                                if path_res.qc.qc_status in ("PASS", "WARNING"):
                                    pathway_status = "PATHWAY_ANALYSIS_SUCCESS"
                                    pathway_metrics = {
                                        "status": "GSEA_SUCCESS",
                                        "contrast_id": primary_contrast.contrast_id,
                                        "qc_status": path_res.qc.qc_status,
                                        "qc_notes": list(path_res.qc.qc_notes),
                                        "ranked_genes_count": ranked_count,
                                        "total_eligible_genes": path_res.ranked_genes.total_eligible_genes,
                                        "go_results_count": go_count,
                                        "kegg_results_count": kegg_count,
                                        "go_significant_pathways_count": path_res.qc.go_significant_pathways_count,
                                        "kegg_significant_pathways_count": path_res.qc.kegg_significant_pathways_count,
                                        "significant_pathways": sig_pathways,
                                    }
                                    final_status = "DOWNSTREAM_COMPLETED"
                                else:
                                    pathway_status = f"PATHWAY_ANALYSIS_QC_{path_res.qc.qc_status}"
                                    pathway_metrics = {
                                        "status": f"GSEA_{path_res.qc.qc_status}",
                                        "contrast_id": primary_contrast.contrast_id,
                                        "qc_status": path_res.qc.qc_status,
                                        "qc_notes": list(path_res.qc.qc_notes),
                                        "ranked_genes_count": ranked_count,
                                        "total_eligible_genes": path_res.ranked_genes.total_eligible_genes,
                                        "go_results_count": go_count,
                                        "kegg_results_count": kegg_count,
                                        "significant_pathways": sig_pathways,
                                    }
                                    final_status = "DOWNSTREAM_PARTIAL_PATHWAY_FAILED"
                            except Exception as e:
                                pathway_status = f"PATHWAY_FAILED: {e}"
                                pathway_metrics = {"reason": str(e), "error": str(e)}
                                final_status = "DOWNSTREAM_PARTIAL_PATHWAY_FAILED"
                        else:
                            deseq2_status = f"DESEQ2_{de_res.status.value}"
                            deseq2_metrics = {
                                "status": deseq2_status,
                                "contrast_id": primary_contrast.contrast_id,
                                "reasons": [f.message for f in de_res.findings],
                            }
                            pathway_status = "PATHWAY_BLOCKED_BY_UPSTREAM_DE"
                            pathway_metrics = {"reason": f"DESeq2 status is {deseq2_status}"}
                            final_status = "DOWNSTREAM_PARTIAL_DESEQ2_FAILED"
                    except Exception as e:
                        deseq2_status = f"DESEQ2_ERROR: {e}"
                        deseq2_metrics = {"error": str(e)}
                        pathway_status = "NOT_EVALUATED"
                        final_status = "DOWNSTREAM_FAILED"
            else:
                deseq2_status = "BLOCKED_MISSING_UPSTREAM_DATASETS"
                pathway_status = "NOT_EVALUATED"
                final_status = "DOWNSTREAM_BLOCKED_MISSING_UPSTREAM"

    return DownstreamBridgeResult(
        status=final_status,
        qc_passed=qc_passed,
        qc_metrics=qc_metrics,
        harmonization_metrics=harm_metrics,
        normalization_metrics=norm_metrics,
        deseq2_status=deseq2_status,
        deseq2_metrics=deseq2_metrics,
        pathway_status=pathway_status,
        pathway_metrics=pathway_metrics,
        biological_replicate_warning=replicate_warning,
        metadata_validation_status=metadata_val_status,
        metadata_validation_metrics=metadata_val_metrics,
    )
