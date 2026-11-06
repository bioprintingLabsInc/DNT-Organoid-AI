"""Stage 8D: Generic Execution Entry Point for Chemical–Gene–Pathway Evidence Relationships.

Provides universal, dataset-independent pipelines to:
1. Ingest Stage 8A.1 registries (in-memory or from CSV).
2. Ingest Stage 8C pathway enrichment results (in-memory or from CSVs).
3. Synthesize structured Chemical–Gene–Pathway evidence records.
4. Execute comprehensive QC and invariant auditing.
5. Serialize to governed canonical formats.
"""

from pathlib import Path
from typing import Any, Sequence
import pandas as pd

from src.features.compound_pathway.builder import (
    build_chemical_gene_records,
    build_chemical_pathway_records,
    link_chemical_gene_pathway_evidence,
)
from src.features.compound_pathway.models import (
    ChemicalExposureContext,
    ChemicalGenePathwayEvidenceDataset,
    PathwayEvidenceStatus,
)
from src.features.compound_pathway.qc import (
    EvidenceQCSummary,
    audit_evidence_dataset_qc,
)
from src.features.compound_pathway.serialization import (
    export_evidence_to_csv,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)
from src.features.pathway.models import PathwayEnrichmentResult


def parse_leading_edge_tokens(val: Any) -> tuple[str, ...]:
    """Parse delimited tokens (semicolon, comma, slash) into a clean tuple."""
    if val is None or pd.isna(val):
        return ()
    s = str(val).strip()
    if not s:
        return ()
    # Semicolon or comma or slash delimited
    delim = ";" if ";" in s else ("," if "," in s else "/")
    return tuple(t.strip() for t in s.split(delim) if t.strip())


def load_pathway_results_from_csv(csv_path: str | Path) -> list[PathwayEnrichmentResult]:
    """Load Stage 8C pathway enrichment results from a canonical CSV export."""
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Pathway results CSV not found at '{path}'.")

    df = pd.read_csv(path)
    results: list[PathwayEnrichmentResult] = []

    for _, row in df.iterrows():
        le_ids = parse_leading_edge_tokens(row.get("leading_edge_gene_ids"))
        le_syms = parse_leading_edge_tokens(row.get("leading_edge_symbols"))

        res = PathwayEnrichmentResult(
            contrast_id=str(row["contrast_id"]),
            source_database=str(row["source_database"]),
            pathway_id=str(row["pathway_id"]),
            pathway_name=str(row["pathway_name"]),
            enrichment_score=float(row["enrichment_score"]),
            normalized_enrichment_score=float(row["normalized_enrichment_score"]),
            p_value=float(row["p_value"]),
            adjusted_p_value=float(row["adjusted_p_value"]),
            pathway_size=int(row["pathway_size"]),
            direction=str(row["direction"]),
            leading_edge_gene_ids=le_ids,
            leading_edge_symbols=le_syms,
            ranking_metric=str(row.get("ranking_metric") or "wald_statistic"),
            analysis_version=str(row.get("analysis_version") or "1.0.0"),
            is_significant=bool(row.get("is_significant", row["adjusted_p_value"] < 0.05)),
        )
        results.append(res)

    return results


def load_registry_from_csv(
    csv_path: str | Path,
    registry_id: str | None = None,
    default_treatment_condition: str = "treatment",
    default_control_condition: str = "control",
) -> ConditionGeneResponseRegistry:
    """Load a ConditionGeneResponseRegistry from a canonical Stage 8A CSV export."""
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Registry CSV not found at '{path}'.")

    df = pd.read_csv(path)
    records: list[ConditionGeneResponseRecord] = []

    for _, row in df.iterrows():
        trt_cond = str(row.get("treatment_condition_id") or default_treatment_condition)
        ctrl_raw = str(row.get("matched_control_condition_ids") or default_control_condition)
        ctrl_tuple = tuple(c.strip() for c in ctrl_raw.split(",") if c.strip()) if ctrl_raw else (default_control_condition,)

        rec = ConditionGeneResponseRecord(
            canonical_ensembl_gene_id=str(row["canonical_ensembl_gene_id"]),
            gene_symbol=str(row["gene_symbol"]) if pd.notna(row.get("gene_symbol")) else None,
            study_id=str(row["study_id"]),
            study_family_id=str(row.get("study_family_id") or row["study_id"]),
            independent_study_id=str(row.get("independent_study_id") or row["study_id"]),
            experiment_id=str(row.get("experiment_id") or "exp_default"),
            contrast_id=str(row["contrast_id"]),
            treatment_condition_id=trt_cond,
            matched_control_condition_ids=ctrl_tuple,
            compound_id=str(row["compound_id"]),
            compound_name=str(row["compound_name"]),
            dose_value=float(row["dose_value"]) if pd.notna(row.get("dose_value")) else None,
            dose_unit=str(row["dose_unit"]) if pd.notna(row.get("dose_unit")) else None,
            dose_original=str(row["dose_original"]) if pd.notna(row.get("dose_original")) else None,
            exposure_duration=float(row["exposure_duration"]) if pd.notna(row.get("exposure_duration")) else None,
            exposure_duration_unit=str(row.get("exposure_duration_unit") or "h"),
            developmental_stage=str(row["developmental_stage"]) if pd.notna(row.get("developmental_stage")) else None,
            model_type=str(row["model_type"]) if pd.notna(row.get("model_type")) else None,
            biological_source_ids=(str(row.get("biological_source_ids") or "source_default"),),
            organoid_context_ids=(str(row.get("organoid_context_ids") or "context_default"),),
            treatment_replicate_count=int(row.get("treatment_replicate_count") or 3),
            control_replicate_count=int(row.get("control_replicate_count") or 3),
            treatment_sample_ids=("sample_trt_1",),
            control_sample_ids=("sample_ctrl_1",),
            treatment_biological_replicate_ids=("rep_trt_1",),
            control_biological_replicate_ids=("rep_ctrl_1",),
            base_mean=float(row["base_mean"]) if pd.notna(row.get("base_mean")) else 0.0,
            log2_fold_change=float(row["log2_fold_change"]) if pd.notna(row.get("log2_fold_change")) else None,
            log2_fold_change_se=float(row["log2_fold_change_se"]) if pd.notna(row.get("log2_fold_change_se")) else None,
            wald_statistic=float(row["wald_statistic"]) if pd.notna(row.get("wald_statistic")) else None,
            p_value=float(row["p_value"]) if pd.notna(row.get("p_value")) else None,
            adjusted_p_value=float(row["adjusted_p_value"]) if pd.notna(row.get("adjusted_p_value")) else None,
            response_direction=ResponseDirection(str(row["response_direction"])),
            estimate_status=EstimateStatus(str(row["estimate_status"])),
            missing_reason=str(row["missing_reason"]) if pd.notna(row.get("missing_reason")) else None,
            cross_study_status=CrossStudyStatus(str(row["cross_study_status"])),
            source_contrast_id=str(row["contrast_id"]),
            de_version=str(row.get("de_version") or "1.0.0"),
            registry_version=str(row.get("registry_version") or "1.0.0"),
        )
        records.append(rec)

    reg_id = registry_id or f"reg_loaded_{path.stem}"
    return ConditionGeneResponseRegistry(records=tuple(records), registry_id=reg_id)


def build_stage8d_evidence_pipeline(
    registry: ConditionGeneResponseRegistry,
    pathway_results: Sequence[PathwayEnrichmentResult] | None = None,
    pathway_status: PathwayEvidenceStatus | None = None,
    contrast_id: str | None = None,
    dataset_id: str | None = None,
    output_csv: str | Path | None = None,
) -> tuple[ChemicalGenePathwayEvidenceDataset, EvidenceQCSummary]:
    """Generic, dataset-independent pipeline to generate Stage 8D evidence records.

    Args:
        registry: Validated Stage 8A.1 ConditionGeneResponseRegistry.
        pathway_results: Optional sequence of Stage 8C PathwayEnrichmentResult objects.
        pathway_status: Explicit qualification status if pathway results are absent/ineligible.
        contrast_id: Optional contrast filter.
        dataset_id: Optional dataset identifier.
        output_csv: Optional destination file path for CSV export.

    Returns:
        tuple of (ChemicalGenePathwayEvidenceDataset, EvidenceQCSummary)
    """
    # 1. Build Chemical-Gene records
    cg_records = build_chemical_gene_records(registry, contrast_id=contrast_id)
    if not cg_records:
        raise ValueError(f"No chemical-gene records found in registry for contrast '{contrast_id}'.")

    # 2. Build Chemical-Pathway records (if available)
    first_exp = cg_records[0].exposure
    cp_records: tuple = ()
    if pathway_results:
        # Filter pathway results matching this contrast if contrast_id is specified
        filtered_pw = [pw for pw in pathway_results if (not contrast_id or pw.contrast_id == contrast_id)]
        if filtered_pw:
            cp_records = build_chemical_pathway_records(filtered_pw, first_exp)

    # 3. Link Evidence
    dataset = link_chemical_gene_pathway_evidence(
        chemical_gene_records=cg_records,
        chemical_pathway_records=cp_records,
        pathway_status=pathway_status,
        dataset_id=dataset_id,
    )

    # 4. Audit Quality Control
    qc = audit_evidence_dataset_qc(dataset)

    # 5. Export if requested
    if output_csv:
        export_evidence_to_csv(dataset, output_csv)

    return dataset, qc
