"""Stage 8D: Generic Execution Entry Point for Chemical–Gene–Pathway Evidence Relationships.

Provides universal, dataset-independent pipelines to:
1. Ingest Stage 8A registries losslessly without fabricating experimental metadata.
2. Ingest Stage 8C pathway enrichment results with robust boolean parsing.
3. Synthesize structured Chemical–Gene–Pathway evidence records enforcing contrast integrity.
4. Execute comprehensive QC and invariant auditing per study, contrast, and gene.
5. Serialize to governed canonical formats.
"""

import json
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


def parse_bool_safely(val: Any, default: bool = False) -> bool:
    """Robustly parse boolean values without string-truthiness traps."""
    if val is None or pd.isna(val):
        return default
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    return default


def parse_leading_edge_tokens(val: Any) -> tuple[str, ...]:
    """Parse delimited tokens (semicolon, comma, slash) into a clean tuple."""
    if val is None or pd.isna(val):
        return ()
    s = str(val).strip()
    if not s:
        return ()
    delim = ";" if ";" in s else ("," if "," in s else "/")
    return tuple(t.strip() for t in s.split(delim) if t.strip())


def load_stage8c_pathway_database_version(provenance_path: str | Path) -> str | None:
    """Extract verified pathway database version from Stage 8C provenance JSON."""
    path = Path(provenance_path)
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        ver = data.get("pathway_database_version")
        return str(ver).strip() if ver and str(ver).strip() else None
    except Exception:
        return None


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

        adj_p = float(row["adjusted_p_value"])
        raw_sig = row.get("is_significant")
        if raw_sig is not None and not pd.isna(raw_sig):
            sig = parse_bool_safely(raw_sig, default=(adj_p < 0.05))
        else:
            sig = bool(adj_p < 0.05)

        raw_ver = str(row["analysis_version"]).strip() if pd.notna(row.get("analysis_version")) and str(row.get("analysis_version")).strip() else None
        res = PathwayEnrichmentResult(
            contrast_id=str(row["contrast_id"]),
            source_database=str(row["source_database"]),
            pathway_id=str(row["pathway_id"]),
            pathway_name=str(row["pathway_name"]),
            enrichment_score=float(row["enrichment_score"]),
            normalized_enrichment_score=float(row["normalized_enrichment_score"]),
            p_value=float(row["p_value"]),
            adjusted_p_value=adj_p,
            pathway_size=int(row["pathway_size"]),
            direction=str(row["direction"]),
            leading_edge_gene_ids=le_ids,
            leading_edge_symbols=le_syms,
            ranking_metric=str(row.get("ranking_metric") or "wald_statistic"),
            analysis_version=raw_ver or "clusterProfiler",
            is_significant=sig,
        )
        results.append(res)

    return results


def load_registry_from_csv(
    csv_path: str | Path,
    registry_id: str | None = None,
) -> ConditionGeneResponseRegistry:
    """Load a ConditionGeneResponseRegistry losslessly without fabricating missing metadata.

    Preserves exact upstream values. Missing metadata remains None or empty tuples;
    never assumes sample IDs, biological replicates, or default replicate counts.
    Never fabricates missing experimental identifiers, duration/dose units, or expression measurements.
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Registry CSV not found at '{path}'.")

    df = pd.read_csv(path)
    records: list[ConditionGeneResponseRecord] = []

    for _, row in df.iterrows():
        # Validate required identifiers - explicitly reject scientifically insufficient inputs
        if pd.isna(row.get("study_id")) or not str(row.get("study_id")).strip():
            raise ValueError("Scientifically insufficient input: missing or empty study_id.")
        if pd.isna(row.get("canonical_ensembl_gene_id")) or not str(row.get("canonical_ensembl_gene_id")).strip():
            raise ValueError("Scientifically insufficient input: missing or empty canonical_ensembl_gene_id.")
        if pd.isna(row.get("contrast_id")) or not str(row.get("contrast_id")).strip():
            raise ValueError("Scientifically insufficient input: missing or empty contrast_id.")

        # Contrast & Condition IDs
        cid = str(row["contrast_id"]).strip()
        trt_cond = str(row["treatment_condition_id"]).strip() if pd.notna(row.get("treatment_condition_id")) and str(row.get("treatment_condition_id")).strip() else None
        ctrl_raw = str(row["matched_control_condition_ids"]) if pd.notna(row.get("matched_control_condition_ids")) else ""
        ctrl_tuple = tuple(c.strip() for c in ctrl_raw.split(";") if c.strip()) if ctrl_raw else ()

        # Replicates and Samples - strictly lossless, no invented identities
        trt_samples_raw = str(row["treatment_sample_ids"]) if pd.notna(row.get("treatment_sample_ids")) else ""
        trt_samples = tuple(s.strip() for s in trt_samples_raw.split(";") if s.strip()) if trt_samples_raw else ()

        ctrl_samples_raw = str(row["control_sample_ids"]) if pd.notna(row.get("control_sample_ids")) else ""
        ctrl_samples = tuple(s.strip() for s in ctrl_samples_raw.split(";") if s.strip()) if ctrl_samples_raw else ()

        trt_bio_raw = str(row["treatment_biological_replicate_ids"]) if pd.notna(row.get("treatment_biological_replicate_ids")) else ""
        trt_bio = tuple(b.strip() for b in trt_bio_raw.split(";") if b.strip()) if trt_bio_raw else ()

        ctrl_bio_raw = str(row["control_biological_replicate_ids"]) if pd.notna(row.get("control_biological_replicate_ids")) else ""
        ctrl_bio = tuple(b.strip() for b in ctrl_bio_raw.split(";") if b.strip()) if ctrl_bio_raw else ()

        # Biological replicate count - NEVER infer from sample counts!
        # Missing replicate counts must remain unknown (None) unless independently supported by verified biological replicate IDs
        trt_count = int(row["treatment_replicate_count"]) if pd.notna(row.get("treatment_replicate_count")) else (len(trt_bio) if trt_bio else None)
        ctrl_count = int(row["control_replicate_count"]) if pd.notna(row.get("control_replicate_count")) else (len(ctrl_bio) if ctrl_bio else None)

        bio_src_raw = str(row["biological_source_ids"]) if pd.notna(row.get("biological_source_ids")) else ""
        bio_src = tuple(s.strip() for s in bio_src_raw.split(";") if s.strip()) if bio_src_raw else ()

        ctx_raw = str(row["organoid_context_ids"]) if pd.notna(row.get("organoid_context_ids")) else ""
        ctx_ids = tuple(c.strip() for c in ctx_raw.split(";") if c.strip()) if ctx_raw else ()

        # Missing experimental identifiers must remain None, not fabricated defaults like "UNKNOWN_EXPERIMENT"
        exp_id = str(row["experiment_id"]).strip() if pd.notna(row.get("experiment_id")) and str(row.get("experiment_id")).strip() else None

        # Study family and independent study IDs
        study_family_id = str(row["study_family_id"]).strip() if pd.notna(row.get("study_family_id")) and str(row.get("study_family_id")).strip() else str(row["study_id"]).strip()
        ind_study_id = str(row["independent_study_id"]).strip() if pd.notna(row.get("independent_study_id")) and str(row.get("independent_study_id")).strip() else str(row["study_id"]).strip()

        # Units and measurements - do not default duration unit to "h" or base_mean to 0.0
        exp_dur_unit = str(row["exposure_duration_unit"]).strip() if pd.notna(row.get("exposure_duration_unit")) and str(row.get("exposure_duration_unit")).strip() else None
        dose_unit = str(row["dose_unit"]).strip() if pd.notna(row.get("dose_unit")) and str(row.get("dose_unit")).strip() else None
        base_mean = float(row["base_mean"]) if pd.notna(row.get("base_mean")) else None

        rec = ConditionGeneResponseRecord(
            canonical_ensembl_gene_id=str(row["canonical_ensembl_gene_id"]).strip(),
            gene_symbol=str(row["gene_symbol"]).strip() if pd.notna(row.get("gene_symbol")) and str(row.get("gene_symbol")).strip() else None,
            study_id=str(row["study_id"]).strip(),
            study_family_id=study_family_id,
            independent_study_id=ind_study_id,
            experiment_id=exp_id,
            contrast_id=cid,
            treatment_condition_id=trt_cond or cid,
            matched_control_condition_ids=ctrl_tuple,
            compound_id=str(row["compound_id"]).strip(),
            compound_name=str(row["compound_name"]).strip(),
            dose_value=float(row["dose_value"]) if pd.notna(row.get("dose_value")) else None,
            dose_unit=dose_unit,
            dose_original=str(row["dose_original"]).strip() if pd.notna(row.get("dose_original")) and str(row.get("dose_original")).strip() else None,
            exposure_duration=float(row["exposure_duration"]) if pd.notna(row.get("exposure_duration")) else None,
            exposure_duration_unit=exp_dur_unit,
            developmental_stage=str(row["developmental_stage"]).strip() if pd.notna(row.get("developmental_stage")) and str(row.get("developmental_stage")).strip() else None,
            model_type=str(row["model_type"]).strip() if pd.notna(row.get("model_type")) and str(row.get("model_type")).strip() else None,
            biological_source_ids=bio_src,
            organoid_context_ids=ctx_ids,
            treatment_replicate_count=trt_count,
            control_replicate_count=ctrl_count,
            treatment_sample_ids=trt_samples,
            control_sample_ids=ctrl_samples,
            treatment_biological_replicate_ids=trt_bio,
            control_biological_replicate_ids=ctrl_bio,
            base_mean=base_mean,
            log2_fold_change=float(row["log2_fold_change"]) if pd.notna(row.get("log2_fold_change")) else None,
            log2_fold_change_se=float(row["log2_fold_change_se"]) if pd.notna(row.get("log2_fold_change_se")) else None,
            wald_statistic=float(row["wald_statistic"]) if pd.notna(row.get("wald_statistic")) else None,
            p_value=float(row["p_value"]) if pd.notna(row.get("p_value")) else None,
            adjusted_p_value=float(row["adjusted_p_value"]) if pd.notna(row.get("adjusted_p_value")) else None,
            response_direction=ResponseDirection(str(row["response_direction"]).strip()),
            estimate_status=EstimateStatus(str(row["estimate_status"]).strip()),
            missing_reason=str(row["missing_reason"]).strip() if pd.notna(row.get("missing_reason")) and str(row.get("missing_reason")).strip() else None,
            cross_study_status=CrossStudyStatus(str(row.get("cross_study_status") or CrossStudyStatus.SINGLE_STUDY_ONLY.value).strip()),
            source_contrast_id=cid,
            de_version=str(row["de_version"]).strip() if pd.notna(row.get("de_version")) and str(row.get("de_version")).strip() else None,
            registry_version=str(row["registry_version"]).strip() if pd.notna(row.get("registry_version")) and str(row.get("registry_version")).strip() else None,
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
    pathway_gene_members: dict[str, set[str]] | None = None,
    pathway_database_version: str | None = None,
    provenance_path: str | Path | None = None,
    has_complete_reference: bool = False,
) -> tuple[ChemicalGenePathwayEvidenceDataset, EvidenceQCSummary]:
    """Generic, dataset-independent pipeline to generate Stage 8D evidence records.

    Args:
        registry: Validated Stage 8A ConditionGeneResponseRegistry.
        pathway_results: Optional sequence of Stage 8C PathwayEnrichmentResult objects.
        pathway_status: Explicit qualification status if pathway results are absent/ineligible.
        contrast_id: Optional contrast filter (if None, processes all contrasts in registry).
        dataset_id: Optional dataset identifier.
        output_csv: Optional destination file path for CSV export.
        pathway_gene_members: Optional dictionary mapping pathway_id to member gene IDs/symbols.
        pathway_database_version: Optional explicit pathway database version from verified Stage 8C provenance.
        provenance_path: Optional path to Stage 8C provenance JSON to extract pathway database version.
        has_complete_reference: Whether pathway_gene_members constitutes a complete, exhaustive gene set reference.

    Returns:
        tuple of (ChemicalGenePathwayEvidenceDataset, EvidenceQCSummary)
    """
    # Resolve pathway database version from verified Stage 8C provenance
    pw_db_version = pathway_database_version
    if pw_db_version is None and provenance_path is not None:
        pw_db_version = load_stage8c_pathway_database_version(provenance_path)
    if pw_db_version is None and output_csv is not None:
        out_p = Path(output_csv)
        candidate_stem = out_p.stem.replace("_chemical_gene_pathway_evidence", "")
        cand_path = out_p.parent / f"{candidate_stem}_provenance.json"
        if cand_path.exists():
            pw_db_version = load_stage8c_pathway_database_version(cand_path)

    # 1. Build Chemical-Gene records
    cg_records = build_chemical_gene_records(registry, contrast_id=contrast_id)
    if not cg_records:
        raise ValueError(f"No chemical-gene records found in registry for contrast '{contrast_id}'.")

    # 2. Build Chemical-Pathway records per contrast
    distinct_contrasts = sorted(list({cg.exposure.contrast_id for cg in cg_records}))
    cp_records: list = []

    if pathway_results:
        # Check that every pathway result matches a contrast in the gene records
        for pw in pathway_results:
            if pw.contrast_id not in distinct_contrasts:
                raise ValueError(
                    f"Contrast mismatch: pathway result '{pw.pathway_id}' belongs to contrast "
                    f"'{pw.contrast_id}', which is not present in the evaluated gene records: {distinct_contrasts}."
                )

        # Build pathway records for each contrast
        for c_id in distinct_contrasts:
            matching_genes = [cg for cg in cg_records if cg.exposure.contrast_id == c_id]
            matching_pathways = [pw for pw in pathway_results if pw.contrast_id == c_id]
            if matching_genes and matching_pathways:
                first_exp = matching_genes[0].exposure
                c_pw_recs = build_chemical_pathway_records(
                    matching_pathways,
                    first_exp,
                    pathway_gene_members=pathway_gene_members,
                    pathway_database_version=pw_db_version,
                    has_complete_membership=has_complete_reference,
                )
                cp_records.extend(c_pw_recs)

    # 3. Link Evidence
    dataset = link_chemical_gene_pathway_evidence(
        chemical_gene_records=cg_records,
        chemical_pathway_records=tuple(cp_records),
        pathway_status=pathway_status,
        dataset_id=dataset_id,
        pathway_gene_members=pathway_gene_members,
        has_complete_reference=has_complete_reference,
    )

    # 4. Audit Quality Control per study, contrast, and gene
    qc = audit_evidence_dataset_qc(dataset)

    # Ensure QC failure prevents a dataset from being marked scientifically valid
    if not qc.is_valid:
        notes_str = "; ".join(qc.notes)
        raise ValueError(
            f"Stage 8D Quality Control verification FAILED for dataset '{dataset.dataset_id}'. "
            f"Reasons: {notes_str}"
        )

    # 5. Export if requested and valid
    if output_csv:
        export_evidence_to_csv(dataset, output_csv)

    return dataset, qc
