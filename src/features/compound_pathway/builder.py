"""Stage 8D: Builder functions for Chemical–Gene–Pathway Evidence Relationships.

Implements deterministic, metadata-driven synthesis of:
    Chemical Exposure -> Observed Gene Perturbation -> Pathway Context
preserving all upstream statistics, enforcing contrast-specific integrity,
and distinguishing the four biological gene-pathway membership states.
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


def _derive_evidence_id(study_id: str, contrast_id: str, gene_id: str, pathway_id: str | None) -> str:
    """Generate deterministic 16-character hex evidence identifier with ev_ prefix.

    Derived from SHA-256(study_id::contrast_id::gene_id::pathway_id).
    Guarantees unique IDs across different contrasts sharing the same gene.
    """
    pw_token = pathway_id.strip() if pathway_id and pathway_id.strip() else "NONE"
    raw_key = f"{study_id.strip()}::{contrast_id.strip()}::{gene_id.strip()}::{pw_token}"
    digest = hashlib.sha256(raw_key.encode("utf-8")).hexdigest()[:16]
    return f"ev_{digest}"


def build_chemical_exposure_context(record: ConditionGeneResponseRecord) -> ChemicalExposureContext:
    """Construct an immutable ChemicalExposureContext from a ConditionGeneResponseRecord losslessly."""
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
        treatment_biological_replicate_ids=record.treatment_biological_replicate_ids,
        control_biological_replicate_ids=record.control_biological_replicate_ids,
        treatment_replicate_count=record.treatment_replicate_count,
        control_replicate_count=record.control_replicate_count,
    )


def build_chemical_gene_records(
    registry: ConditionGeneResponseRegistry,
    contrast_id: str | None = None,
) -> tuple[ChemicalGeneRecord, ...]:
    """Convert Stage 8A.1/8A.2 registry records into typed ChemicalGeneRecords.

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
            causality_claim="OBSERVED_TRANSCRIPTOMIC_RESPONSE",
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
    pathway_gene_members: dict[str, set[str]] | None = None,
) -> tuple[ChemicalPathwayRecord, ...]:
    """Convert Stage 8C GSEA results into typed ChemicalPathwayRecords.

    Enforces contrast-specific matching: rejects mismatched pathway results.
    """
    cp_records: list[ChemicalPathwayRecord] = []
    for pw in pathway_results:
        # Enforce contrast-specific identity match
        if pw.contrast_id != exposure.contrast_id:
            raise ValueError(
                f"Contrast identity mismatch: pathway result '{pw.pathway_id}' belongs to "
                f"contrast '{pw.contrast_id}', but exposure context is '{exposure.contrast_id}'."
            )

        # Attach complete reference members if available
        all_ids: tuple[str, ...] = ()
        all_syms: tuple[str, ...] = ()
        if pathway_gene_members and pw.pathway_id in pathway_gene_members:
            members = pathway_gene_members[pw.pathway_id]
            all_ids = tuple(m for m in members if m.startswith("ENSG"))
            all_syms = tuple(m for m in members if not m.startswith("ENSG"))

        rec = ChemicalPathwayRecord(
            exposure=exposure,
            pathway_id=pw.pathway_id,
            pathway_name=pw.pathway_name,
            source_database=pw.source_database,
            database_version="1.0.0",
            enrichment_score=pw.enrichment_score,
            normalized_enrichment_score=pw.normalized_enrichment_score,
            p_value=pw.p_value,
            adjusted_p_value=pw.adjusted_p_value,
            pathway_size=pw.pathway_size,
            pathway_direction=pw.direction,
            is_significant=pw.is_significant,
            leading_edge_gene_ids=pw.leading_edge_gene_ids,
            leading_edge_symbols=pw.leading_edge_symbols,
            all_member_gene_ids=all_ids,
            all_member_symbols=all_syms,
            ranking_metric=pw.ranking_metric,
            analysis_version=pw.analysis_version,
            scientific_limitation=(
                "Statistically evaluated pathway enrichment from preranked GSEA. "
                "Reflects coordinated transcriptional perturbation across the gene set."
            ),
        )
        cp_records.append(rec)
    return tuple(cp_records)


def determine_gene_pathway_membership(
    gene_id: str,
    gene_symbol: str | None,
    pathway: ChemicalPathwayRecord,
    pathway_gene_members: dict[str, set[str]] | None = None,
) -> tuple[GenePathwayMembershipType, bool]:
    """Determine the four-state biological gene-pathway participation.

    Distinguishes:
    1. LEADING_EDGE_MEMBER: Verified core driving member in GSEA.
    2. NON_LEADING_EDGE_MEMBER: Verified member of pathway gene set outside the leading edge.
    3. CONFIRMED_NON_MEMBER: Complete versioned gene set reference confirms gene is NOT a member.
    4. UNKNOWN_OR_UNAVAILABLE: Complete gene set reference was unavailable or incomplete.

    Returns:
        tuple of (GenePathwayMembershipType, is_leading_edge: bool)
    """
    le_tokens = {g.strip().upper() for g in pathway.leading_edge_gene_ids + pathway.leading_edge_symbols if g}
    gid_clean = gene_id.strip().upper() if gene_id else ""
    sym_clean = gene_symbol.strip().upper() if gene_symbol else ""

    # 1. Check verified leading edge membership
    if (gid_clean and gid_clean in le_tokens) or (sym_clean and sym_clean in le_tokens):
        return GenePathwayMembershipType.LEADING_EDGE_MEMBER, True

    # 2. Check complete reference membership if available
    ref_tokens: set[str] = set()
    if pathway.all_member_gene_ids or pathway.all_member_symbols:
        ref_tokens = {g.strip().upper() for g in pathway.all_member_gene_ids + pathway.all_member_symbols if g}
    elif pathway_gene_members and pathway.pathway_id in pathway_gene_members:
        ref_tokens = {g.strip().upper() for g in pathway_gene_members[pathway.pathway_id] if g}

    if ref_tokens:
        if (gid_clean and gid_clean in ref_tokens) or (sym_clean and sym_clean in ref_tokens):
            return GenePathwayMembershipType.NON_LEADING_EDGE_MEMBER, False
        else:
            return GenePathwayMembershipType.CONFIRMED_NON_MEMBER, False

    # 3. Reference is unavailable or incomplete
    return GenePathwayMembershipType.UNKNOWN_OR_UNAVAILABLE, False


def link_chemical_gene_pathway_evidence(
    chemical_gene_records: Sequence[ChemicalGeneRecord],
    chemical_pathway_records: Sequence[ChemicalPathwayRecord] | None = None,
    pathway_status: PathwayEvidenceStatus | None = None,
    dataset_id: str | None = None,
    pathway_gene_members: dict[str, set[str]] | None = None,
) -> ChemicalGenePathwayEvidenceDataset:
    """Synthesize structured Chemical–Gene–Pathway evidence records with strict contrast enforcement.

    Enforces:
    1. Law of Conservation of the Gene Universe: 100% of input genes are preserved.
    2. Strict Contrast Integrity: Rejects mismatched pathway results across contrasts.
    3. Four-State Membership: Distinguishes leading-edge, non-leading-edge member, confirmed non-member,
       and unavailable membership.
    4. Decoupling: Separates gene fold-change sign from pathway direction and pathway significance.
    """
    if not chemical_gene_records:
        raise ValueError("Cannot link evidence with empty chemical_gene_records.")

    pw_records_seq = tuple(chemical_pathway_records or ())

    # Group gene records by (study_id, contrast_id)
    genes_by_contrast: dict[tuple[str, str], list[ChemicalGeneRecord]] = {}
    for cg in chemical_gene_records:
        key = (cg.exposure.study_id, cg.exposure.contrast_id)
        genes_by_contrast.setdefault(key, []).append(cg)

    # Group pathway records by (study_id, contrast_id)
    pathways_by_contrast: dict[tuple[str, str], list[ChemicalPathwayRecord]] = {}
    for pw in pw_records_seq:
        key = (pw.exposure.study_id, pw.exposure.contrast_id)
        # Check if pathway's contrast exists in the input gene records
        if key not in genes_by_contrast:
            raise ValueError(
                f"Mismatched contrast rejection: Pathway result '{pw.pathway_id}' belongs to contrast "
                f"'{pw.exposure.contrast_id}' in study '{pw.exposure.study_id}', which does not exist "
                f"in the provided gene records."
            )
        pathways_by_contrast.setdefault(key, []).append(pw)

    evidence_records: list[ChemicalGenePathwayEvidenceRecord] = []

    # Process per (study_id, contrast_id)
    for (study_id, contrast_id), contrast_genes in genes_by_contrast.items():
        contrast_pathways = pathways_by_contrast.get((study_id, contrast_id), [])
        has_pathways = len(contrast_pathways) > 0

        # Determine qualification status for this contrast
        if pathway_status is not None:
            current_status = pathway_status
        elif has_pathways:
            current_status = PathwayEvidenceStatus.AVAILABLE
        else:
            current_status = PathwayEvidenceStatus.NOT_AVAILABLE

        # Build high-performance inverted index for fast gene-to-pathway mapping
        gene_to_leading_edge: dict[str, list[ChemicalPathwayRecord]] = {}
        gene_to_outside_member: dict[str, list[ChemicalPathwayRecord]] = {}

        if current_status == PathwayEvidenceStatus.AVAILABLE and has_pathways:
            for pw in contrast_pathways:
                # 1. Index leading-edge genes
                le_tokens = {g.strip().upper() for g in pw.leading_edge_gene_ids + pw.leading_edge_symbols if g}
                for tok in le_tokens:
                    gene_to_leading_edge.setdefault(tok, []).append(pw)

                # 2. Index outside members if complete reference is available
                ref_tokens: set[str] = set()
                if pw.all_member_gene_ids or pw.all_member_symbols:
                    ref_tokens = {g.strip().upper() for g in pw.all_member_gene_ids + pw.all_member_symbols if g}
                elif pathway_gene_members and pw.pathway_id in pathway_gene_members:
                    ref_tokens = {g.strip().upper() for g in pathway_gene_members[pw.pathway_id] if g}

                outside_tokens = ref_tokens - le_tokens
                for tok in outside_tokens:
                    gene_to_outside_member.setdefault(tok, []).append(pw)

        for cg in contrast_genes:
            exp = cg.exposure
            gid_clean = cg.canonical_ensembl_gene_id.strip().upper()
            sym_clean = cg.gene_symbol.strip().upper() if cg.gene_symbol else None

            matched_records: list[ChemicalGenePathwayEvidenceRecord] = []

            if current_status == PathwayEvidenceStatus.AVAILABLE and has_pathways:
                # 1. Collect leading-edge matches
                le_pws = list(gene_to_leading_edge.get(gid_clean, []))
                if sym_clean and sym_clean in gene_to_leading_edge:
                    for p in gene_to_leading_edge[sym_clean]:
                        if p not in le_pws:
                            le_pws.append(p)

                for pw in le_pws:
                    ev_id = _derive_evidence_id(study_id, contrast_id, cg.canonical_ensembl_gene_id, pw.pathway_id)
                    rec = ChemicalGenePathwayEvidenceRecord(
                        evidence_id=ev_id,
                        study_id=exp.study_id,
                        contrast_id=exp.contrast_id,
                        independent_study_id=exp.independent_study_id,
                        compound_id=exp.compound_id,
                        compound_name=exp.compound_name,
                        treatment_condition_id=exp.treatment_condition_id,
                        matched_control_condition_ids=exp.matched_control_condition_ids,
                        dose_value=exp.dose_value,
                        dose_unit=exp.dose_unit,
                        exposure_duration=exp.exposure_duration,
                        exposure_duration_unit=exp.exposure_duration_unit,
                        developmental_stage=exp.developmental_stage,
                        model_type=exp.model_type,
                        canonical_ensembl_gene_id=cg.canonical_ensembl_gene_id,
                        gene_symbol=cg.gene_symbol,
                        base_mean=cg.base_mean,
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
                        enrichment_score=pw.enrichment_score,
                        pathway_normalized_enrichment_score=pw.normalized_enrichment_score,
                        pathway_p_value=pw.p_value,
                        pathway_adjusted_p_value=pw.adjusted_p_value,
                        pathway_size=pw.pathway_size,
                        pathway_direction=pw.pathway_direction,
                        is_pathway_significant=pw.is_significant,
                        gene_membership_type=GenePathwayMembershipType.LEADING_EDGE_MEMBER,
                        is_leading_edge=True,
                        evidence_type=EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION,
                        is_direct_target_interaction=False,
                        causality_claim="OBSERVED_TRANSCRIPTOMIC_RESPONSE",
                        scientific_limitations=(
                            "Gene observed as leading-edge driver in statistically evaluated pathway. "
                            "Does NOT imply direct chemical-gene binding or directional pathway causality."
                        ),
                        provenance_notes=f"Pathway source: {pw.source_database} ({pw.analysis_version})",
                    )
                    matched_records.append(rec)

                # 2. Collect outside member matches
                outside_pws = list(gene_to_outside_member.get(gid_clean, []))
                if sym_clean and sym_clean in gene_to_outside_member:
                    for p in gene_to_outside_member[sym_clean]:
                        if p not in outside_pws and p not in le_pws:
                            outside_pws.append(p)

                for pw in outside_pws:
                    ev_id = _derive_evidence_id(study_id, contrast_id, cg.canonical_ensembl_gene_id, pw.pathway_id)
                    rec = ChemicalGenePathwayEvidenceRecord(
                        evidence_id=ev_id,
                        study_id=exp.study_id,
                        contrast_id=exp.contrast_id,
                        independent_study_id=exp.independent_study_id,
                        compound_id=exp.compound_id,
                        compound_name=exp.compound_name,
                        treatment_condition_id=exp.treatment_condition_id,
                        matched_control_condition_ids=exp.matched_control_condition_ids,
                        dose_value=exp.dose_value,
                        dose_unit=exp.dose_unit,
                        exposure_duration=exp.exposure_duration,
                        exposure_duration_unit=exp.exposure_duration_unit,
                        developmental_stage=exp.developmental_stage,
                        model_type=exp.model_type,
                        canonical_ensembl_gene_id=cg.canonical_ensembl_gene_id,
                        gene_symbol=cg.gene_symbol,
                        base_mean=cg.base_mean,
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
                        enrichment_score=pw.enrichment_score,
                        pathway_normalized_enrichment_score=pw.normalized_enrichment_score,
                        pathway_p_value=pw.p_value,
                        pathway_adjusted_p_value=pw.adjusted_p_value,
                        pathway_size=pw.pathway_size,
                        pathway_direction=pw.pathway_direction,
                        is_pathway_significant=pw.is_significant,
                        gene_membership_type=GenePathwayMembershipType.NON_LEADING_EDGE_MEMBER,
                        is_leading_edge=False,
                        evidence_type=EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION,
                        is_direct_target_interaction=False,
                        causality_claim="OBSERVED_TRANSCRIPTOMIC_RESPONSE",
                        scientific_limitations=(
                            "Gene observed as pathway gene set member outside the leading edge. "
                            "Does NOT imply direct chemical-gene binding or directional pathway causality."
                        ),
                        provenance_notes=f"Pathway source: {pw.source_database} ({pw.analysis_version})",
                    )
                    matched_records.append(rec)

            if matched_records:
                evidence_records.extend(matched_records)
            else:
                # Gene conserved without matching evaluated pathway membership
                ev_id = _derive_evidence_id(study_id, contrast_id, cg.canonical_ensembl_gene_id, None)
                if current_status == PathwayEvidenceStatus.AVAILABLE:
                    pw_stat = PathwayEvidenceStatus.NOT_ANNOTATED
                    mem_type = GenePathwayMembershipType.CONFIRMED_NON_MEMBER if pathway_gene_members else GenePathwayMembershipType.UNKNOWN_OR_UNAVAILABLE
                    limitation_text = "Observed transcriptomic perturbation. Gene is not an annotated member or leading-edge driver of evaluated pathways."
                    prov_note = "Conserved gene universe record; not participating in evaluated pathways"
                else:
                    pw_stat = current_status
                    mem_type = GenePathwayMembershipType.UNKNOWN_OR_UNAVAILABLE
                    limitation_text = "Observed transcriptomic perturbation. Pathway analysis not performed or unavailable for this contrast."
                    prov_note = f"Conserved gene universe record; pathway status: {current_status}"

                rec = ChemicalGenePathwayEvidenceRecord(
                    evidence_id=ev_id,
                    study_id=exp.study_id,
                    contrast_id=exp.contrast_id,
                    independent_study_id=exp.independent_study_id,
                    compound_id=exp.compound_id,
                    compound_name=exp.compound_name,
                    treatment_condition_id=exp.treatment_condition_id,
                    matched_control_condition_ids=exp.matched_control_condition_ids,
                    dose_value=exp.dose_value,
                    dose_unit=exp.dose_unit,
                    exposure_duration=exp.exposure_duration,
                    exposure_duration_unit=exp.exposure_duration_unit,
                    developmental_stage=exp.developmental_stage,
                    model_type=exp.model_type,
                    canonical_ensembl_gene_id=cg.canonical_ensembl_gene_id,
                    gene_symbol=cg.gene_symbol,
                    base_mean=cg.base_mean,
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
                    enrichment_score=None,
                    pathway_normalized_enrichment_score=None,
                    pathway_p_value=None,
                    pathway_adjusted_p_value=None,
                    pathway_size=None,
                    pathway_direction=None,
                    is_pathway_significant=None,
                    gene_membership_type=mem_type,
                    is_leading_edge=False,
                    evidence_type=EvidenceType.INACTIVE_OR_UNANNOTATED_PERTURBATION,
                    is_direct_target_interaction=False,
                    causality_claim="OBSERVED_TRANSCRIPTOMIC_RESPONSE",
                    scientific_limitations=limitation_text,
                    provenance_notes=prov_note,
                )
                evidence_records.append(rec)

    ds_id = dataset_id or f"evds_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
    return ChemicalGenePathwayEvidenceDataset(
        dataset_id=ds_id,
        records=tuple(evidence_records),
        chemical_gene_records=tuple(chemical_gene_records),
        chemical_pathway_records=pw_records_seq,
        created_at=datetime.now(timezone.utc).isoformat(),
        analysis_version="1.0.0",
    )
