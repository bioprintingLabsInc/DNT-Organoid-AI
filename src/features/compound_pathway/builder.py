"""Stage 8D: Builder functions for Chemical–Gene–Pathway Evidence Relationships.

Implements deterministic, metadata-driven synthesis of:
    Chemical Exposure -> Observed Gene Perturbation -> Pathway Context
preserving all upstream statistics and maintaining universal chemical/dataset independence.
"""

from collections.abc import Sequence
from datetime import datetime, timezone
import hashlib
from typing import Any

from src.features.compound_pathway.models import (
    ChemicalExposureContext,
    ChemicalGenePathwayEvidenceDataset,
    ChemicalGenePathwayEvidenceRecord,
    ChemicalGeneRecord,
    ChemicalPathwayRecord,
    EvidenceType,
    GenePathwayMembershipType,
    PathwayEvidenceStatus,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
)
from src.features.pathway.models import PathwayEnrichmentResult


def build_chemical_exposure_context(record: ConditionGeneResponseRecord) -> ChemicalExposureContext:
    """Construct an immutable ChemicalExposureContext from a ConditionGeneResponseRecord."""
    return ChemicalExposureContext(
        compound_id=record.compound_id,
        compound_name=record.compound_name,
        dose_value=record.dose_value,
        dose_unit=record.dose_unit,
        dose_original=record.dose_original,
        exposure_duration=record.exposure_duration,
        exposure_duration_unit=record.exposure_duration_unit,
        developmental_stage=record.developmental_stage,
        model_type=record.model_type,
        biological_source_ids=record.biological_source_ids,
        organoid_context_ids=record.organoid_context_ids,
        study_id=record.study_id,
        study_family_id=record.study_family_id,
        independent_study_id=record.independent_study_id,
        experiment_id=record.experiment_id,
        contrast_id=record.contrast_id,
        treatment_condition_id=record.treatment_condition_id,
        matched_control_condition_ids=record.matched_control_condition_ids,
        treatment_sample_ids=record.treatment_sample_ids,
        control_sample_ids=record.control_sample_ids,
        treatment_replicate_count=record.treatment_replicate_count,
        control_replicate_count=record.control_replicate_count,
    )


def build_chemical_gene_records(
    registry: ConditionGeneResponseRegistry,
    contrast_id: str | None = None,
) -> tuple[ChemicalGeneRecord, ...]:
    """Convert Stage 8A.1 registry records into typed ChemicalGeneRecords.

    Preserves 100% of underlying differential expression statistics without modification.
    """
    records = registry.get_records_for_contrast(contrast_id) if contrast_id else registry.records
    cg_records: list[ChemicalGeneRecord] = []

    for r in records:
        exp = build_chemical_exposure_context(r)
        cg_rec = ChemicalGeneRecord(
            exposure=exp,
            canonical_ensembl_gene_id=r.canonical_ensembl_gene_id,
            gene_symbol=r.gene_symbol,
            entrez_gene_id=None,
            base_mean=r.base_mean,
            log2_fold_change=r.log2_fold_change,
            log2_fold_change_se=r.log2_fold_change_se,
            wald_statistic=r.wald_statistic,
            p_value=r.p_value,
            adjusted_p_value=r.adjusted_p_value,
            response_direction=r.response_direction,
            estimate_status=r.estimate_status,
            evidence_type=EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION,
            is_direct_target_interaction=False,
            scientific_limitation=(
                "Observed transcriptomic perturbation following chemical exposure. "
                "Does NOT constitute verified direct binding or direct biochemical target interaction."
            ),
            upstream_provenance={
                "source_contrast_id": r.source_contrast_id,
                "de_version": r.de_version,
                "registry_version": r.registry_version,
            },
        )
        cg_records.append(cg_rec)

    return tuple(cg_records)


def build_chemical_pathway_records(
    pathway_results: Sequence[PathwayEnrichmentResult],
    exposure: ChemicalExposureContext,
    database_version: str | None = None,
) -> tuple[ChemicalPathwayRecord, ...]:
    """Convert Stage 8C pathway enrichment results into typed ChemicalPathwayRecords."""
    cp_records: list[ChemicalPathwayRecord] = []

    for pw in pathway_results:
        cp_rec = ChemicalPathwayRecord(
            exposure=exposure,
            pathway_id=pw.pathway_id,
            pathway_name=pw.pathway_name,
            source_database=pw.source_database,
            database_version=database_version,
            enrichment_score=pw.enrichment_score,
            normalized_enrichment_score=pw.normalized_enrichment_score,
            p_value=pw.p_value,
            adjusted_p_value=pw.adjusted_p_value,
            pathway_size=pw.pathway_size,
            pathway_direction=pw.direction,
            is_significant=pw.is_significant,
            leading_edge_gene_ids=tuple(pw.leading_edge_gene_ids),
            leading_edge_symbols=tuple(pw.leading_edge_symbols),
            ranking_metric=pw.ranking_metric,
            analysis_version=pw.analysis_version,
            scientific_limitation=(
                "Statistically evaluated pathway enrichment from preranked GSEA. "
                "Reflects coordinated transcriptional modulation across the gene set."
            ),
        )
        cp_records.append(cp_rec)

    return tuple(cp_records)


def _derive_evidence_id(contrast_id: str, gene_id: str, pathway_id: str | None) -> str:
    """Generate deterministic unique identifier for a chemical-gene-pathway evidence relationship."""
    raw = f"{contrast_id}::{gene_id}::{pathway_id or 'NONE'}"
    h = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]
    return f"ev_{h}"


def link_chemical_gene_pathway_evidence(
    chemical_gene_records: Sequence[ChemicalGeneRecord],
    chemical_pathway_records: Sequence[ChemicalPathwayRecord] | None = None,
    pathway_status: PathwayEvidenceStatus | None = None,
    dataset_id: str | None = None,
) -> ChemicalGenePathwayEvidenceDataset:
    """Synthesize tripartite Chemical–Gene–Pathway evidence records.

    Law of Conservation of the Gene Universe:
    Every input gene is preserved in the evidence dataset.
    - If pathway enrichment is available:
        - Genes belonging to pathways (leading edge) receive linked pathway records.
        - Genes not in evaluated pathway leading edges receive unannotated baseline records.
    - If pathway enrichment is unavailable:
        - All genes are preserved with status NOT_AVAILABLE and explicit limitations.
    """
    if not chemical_gene_records:
        raise ValueError("Cannot construct evidence dataset: chemical_gene_records is empty.")

    first_exp = chemical_gene_records[0].exposure
    cid = dataset_id or f"evds_{first_exp.contrast_id}"
    timestamp = datetime.now(timezone.utc).isoformat()

    pw_records = tuple(chemical_pathway_records) if chemical_pathway_records else ()
    status = pathway_status or (PathwayEvidenceStatus.AVAILABLE if pw_records else PathwayEvidenceStatus.NOT_AVAILABLE)

    # Build index of pathways by leading edge genes (supports canonical Ensembl IDs and gene symbols)
    gene_to_pathways: dict[str, list[ChemicalPathwayRecord]] = {}
    if status == PathwayEvidenceStatus.AVAILABLE and pw_records:
        for pw in pw_records:
            # Match by leading edge gene IDs (Ensembl or Entrez)
            for gid in pw.leading_edge_gene_ids:
                if gid:
                    gene_to_pathways.setdefault(gid.strip(), []).append(pw)
            # Match by symbols
            for sym in pw.leading_edge_symbols:
                if sym:
                    gene_to_pathways.setdefault(sym.strip().upper(), []).append(pw)

    evidence_records: list[ChemicalGenePathwayEvidenceRecord] = []

    for cg in chemical_gene_records:
        exp = cg.exposure
        gid = cg.canonical_ensembl_gene_id
        sym = cg.gene_symbol.strip().upper() if cg.gene_symbol else None

        # Find matching pathways for this gene
        matched_pathways: list[ChemicalPathwayRecord] = []
        if status == PathwayEvidenceStatus.AVAILABLE:
            if gid in gene_to_pathways:
                matched_pathways.extend(gene_to_pathways[gid])
            if sym and sym in gene_to_pathways:
                for p in gene_to_pathways[sym]:
                    if p not in matched_pathways:
                        matched_pathways.append(p)

        if matched_pathways:
            for pw in matched_pathways:
                ev_id = _derive_evidence_id(exp.contrast_id, gid, pw.pathway_id)
                rec = ChemicalGenePathwayEvidenceRecord(
                    evidence_id=ev_id,
                    contrast_id=exp.contrast_id,
                    study_id=exp.study_id,
                    independent_study_id=exp.independent_study_id,
                    compound_id=exp.compound_id,
                    compound_name=exp.compound_name,
                    dose_value=exp.dose_value,
                    dose_unit=exp.dose_unit,
                    exposure_duration=exp.exposure_duration,
                    exposure_duration_unit=exp.exposure_duration_unit,
                    developmental_stage=exp.developmental_stage,
                    model_type=exp.model_type,
                    canonical_ensembl_gene_id=gid,
                    gene_symbol=cg.gene_symbol,
                    log2_fold_change=cg.log2_fold_change,
                    log2_fold_change_se=cg.log2_fold_change_se,
                    wald_statistic=cg.wald_statistic,
                    p_value=cg.p_value,
                    adjusted_p_value=cg.adjusted_p_value,
                    gene_response_direction=cg.response_direction,
                    gene_estimate_status=cg.estimate_status,
                    pathway_evidence_status=PathwayEvidenceStatus.AVAILABLE,
                    pathway_id=pw.pathway_id,
                    pathway_name=pw.pathway_name,
                    source_database=pw.source_database,
                    database_version=pw.database_version,
                    pathway_normalized_enrichment_score=pw.normalized_enrichment_score,
                    pathway_p_value=pw.p_value,
                    pathway_adjusted_p_value=pw.adjusted_p_value,
                    pathway_direction=pw.pathway_direction,
                    is_pathway_significant=pw.is_significant,
                    gene_membership_type=GenePathwayMembershipType.LEADING_EDGE_MEMBER,
                    is_leading_edge=True,
                    evidence_type=EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION,
                    is_direct_target_interaction=False,
                    scientific_limitations=(
                        "Gene observed as leading-edge driver in statistically enriched pathway. "
                        "Does NOT imply direct chemical-gene binding or directional pathway causality."
                    ),
                    provenance_notes=f"Pathway source: {pw.source_database} ({pw.analysis_version})",
                )
                evidence_records.append(rec)
        else:
            # Gene conserved without matching pathway leading edge or pathway evidence unavailable
            pw_stat = status
            limitations = (
                "Observed transcriptomic perturbation. Pathway analysis not performed or unavailable for this contrast."
                if status != PathwayEvidenceStatus.AVAILABLE
                else "Observed transcriptomic perturbation. Gene is not a leading-edge driver of enriched pathways."
            )
            ev_id = _derive_evidence_id(exp.contrast_id, gid, None)
            rec = ChemicalGenePathwayEvidenceRecord(
                evidence_id=ev_id,
                contrast_id=exp.contrast_id,
                study_id=exp.study_id,
                independent_study_id=exp.independent_study_id,
                compound_id=exp.compound_id,
                compound_name=exp.compound_name,
                dose_value=exp.dose_value,
                dose_unit=exp.dose_unit,
                exposure_duration=exp.exposure_duration,
                exposure_duration_unit=exp.exposure_duration_unit,
                developmental_stage=exp.developmental_stage,
                model_type=exp.model_type,
                canonical_ensembl_gene_id=gid,
                gene_symbol=cg.gene_symbol,
                log2_fold_change=cg.log2_fold_change,
                log2_fold_change_se=cg.log2_fold_change_se,
                wald_statistic=cg.wald_statistic,
                p_value=cg.p_value,
                adjusted_p_value=cg.adjusted_p_value,
                gene_response_direction=cg.response_direction,
                gene_estimate_status=cg.estimate_status,
                pathway_evidence_status=pw_stat,
                pathway_id=None,
                pathway_name=None,
                source_database=None,
                database_version=None,
                pathway_normalized_enrichment_score=None,
                pathway_p_value=None,
                pathway_adjusted_p_value=None,
                pathway_direction=None,
                is_pathway_significant=None,
                gene_membership_type=GenePathwayMembershipType.UNANNOTATED,
                is_leading_edge=False,
                evidence_type=EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION,
                is_direct_target_interaction=False,
                scientific_limitations=limitations,
                provenance_notes="Conserved gene universe record",
            )
            evidence_records.append(rec)

    return ChemicalGenePathwayEvidenceDataset(
        dataset_id=cid,
        records=tuple(evidence_records),
        chemical_gene_records=tuple(chemical_gene_records),
        chemical_pathway_records=pw_records,
        created_at=timestamp,
        analysis_version="1.0.0",
    )
