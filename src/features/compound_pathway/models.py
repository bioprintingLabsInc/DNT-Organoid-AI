"""Stage 8D: Data Models for Chemical–Gene–Pathway Evidence Relationships.

Defines immutable, typed dataclasses and controlled vocabularies connecting:
    Chemical Exposure -> Observed Gene Perturbation -> Pathway Context
maintaining exact upstream statistics, clear membership boundaries, and strict
safeguards against unverified physical binding or causality claims.
"""

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any

from src.features.cross_study.models import (
    EstimateStatus,
    ResponseDirection,
)


class EvidenceType(StrEnum):
    """Governed taxonomy of experimental evidence relationships."""

    OBSERVED_TRANSCRIPTOMIC_PERTURBATION = "OBSERVED_TRANSCRIPTOMIC_PERTURBATION"
    PATHWAY_COLLECTIVE_PERTURBATION = "PATHWAY_COLLECTIVE_PERTURBATION"
    INACTIVE_OR_UNANNOTATED_PERTURBATION = "INACTIVE_OR_UNANNOTATED_PERTURBATION"


class PathwayEvidenceStatus(StrEnum):
    """Qualification status of biological pathway evidence for a contrast."""

    AVAILABLE = "AVAILABLE"
    NOT_ANNOTATED = "NOT_ANNOTATED"
    NOT_AVAILABLE = "NOT_AVAILABLE"
    EVALUATION_INELIGIBLE = "EVALUATION_INELIGIBLE"
    INSUFFICIENT_GENES = "INSUFFICIENT_GENES"


class GenePathwayMembershipType(StrEnum):
    """Nature of gene participation in a biological pathway.

    Distinguishes:
    1. LEADING_EDGE_MEMBER: Verified core driving member in preranked GSEA.
    2. NON_LEADING_EDGE_MEMBER: Verified member of the pathway gene set outside the leading edge.
    3. CONFIRMED_NON_MEMBER: Verified complete gene set reference confirms gene is NOT a member.
    4. UNKNOWN_OR_UNAVAILABLE: Complete gene set reference was unavailable, or pathway was not evaluated.
    """

    LEADING_EDGE_MEMBER = "LEADING_EDGE_MEMBER"
    NON_LEADING_EDGE_MEMBER = "NON_LEADING_EDGE_MEMBER"
    CONFIRMED_NON_MEMBER = "CONFIRMED_NON_MEMBER"
    UNKNOWN_OR_UNAVAILABLE = "UNKNOWN_OR_UNAVAILABLE"


@dataclass(frozen=True, slots=True)
class ChemicalExposureContext:
    """Comprehensive chemical and experimental exposure metadata for a single contrast."""

    compound_id: str
    compound_name: str
    dose_value: float | None
    dose_unit: str | None
    dose_original: str | None
    exposure_duration: float | None
    exposure_duration_unit: str | None
    developmental_stage: str | None
    model_type: str | None
    biological_source_ids: tuple[str, ...] = ()
    organoid_context_ids: tuple[str, ...] = ()
    study_id: str = ""
    study_family_id: str | None = None
    independent_study_id: str | None = None
    experiment_id: str | None = None
    contrast_id: str = ""
    treatment_condition_id: str | None = None
    matched_control_condition_ids: tuple[str, ...] = ()
    treatment_sample_ids: tuple[str, ...] = ()
    control_sample_ids: tuple[str, ...] = ()
    treatment_biological_replicate_ids: tuple[str, ...] = ()
    control_biological_replicate_ids: tuple[str, ...] = ()
    treatment_replicate_count: int | None = None
    control_replicate_count: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "dose_value": self.dose_value,
            "dose_unit": self.dose_unit,
            "dose_original": self.dose_original,
            "exposure_duration": self.exposure_duration,
            "exposure_duration_unit": self.exposure_duration_unit,
            "developmental_stage": self.developmental_stage,
            "model_type": self.model_type,
            "biological_source_ids": list(self.biological_source_ids),
            "organoid_context_ids": list(self.organoid_context_ids),
            "study_id": self.study_id,
            "study_family_id": self.study_family_id,
            "independent_study_id": self.independent_study_id,
            "experiment_id": self.experiment_id,
            "contrast_id": self.contrast_id,
            "treatment_condition_id": self.treatment_condition_id,
            "matched_control_condition_ids": list(self.matched_control_condition_ids),
            "treatment_sample_ids": list(self.treatment_sample_ids),
            "control_sample_ids": list(self.control_sample_ids),
            "treatment_biological_replicate_ids": list(self.treatment_biological_replicate_ids),
            "control_biological_replicate_ids": list(self.control_biological_replicate_ids),
            "treatment_replicate_count": self.treatment_replicate_count,
            "control_replicate_count": self.control_replicate_count,
        }


@dataclass(frozen=True, slots=True)
class ChemicalGeneRecord:
    """Observed chemical-to-gene transcriptional response relationship."""

    exposure: ChemicalExposureContext
    canonical_ensembl_gene_id: str
    gene_symbol: str | None
    entrez_gene_id: str | None
    base_mean: float | None = None
    log2_fold_change: float | None = None
    log2_fold_change_se: float | None = None
    wald_statistic: float | None = None
    p_value: float | None = None
    adjusted_p_value: float | None = None
    response_direction: ResponseDirection = ResponseDirection.NO_DIRECTION
    estimate_status: EstimateStatus = EstimateStatus.OK
    evidence_type: EvidenceType = EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION
    is_direct_target_interaction: bool = False
    causality_claim: str = "OBSERVED_TRANSCRIPTOMIC_RESPONSE"
    scientific_limitation: str = (
        "Observed transcriptomic perturbation following chemical exposure. "
        "Does NOT constitute verified direct binding or direct biochemical target interaction."
    )
    upstream_provenance: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["exposure"] = self.exposure.to_dict()
        d["response_direction"] = str(self.response_direction)
        d["estimate_status"] = str(self.estimate_status)
        d["evidence_type"] = str(self.evidence_type)
        return d


@dataclass(frozen=True, slots=True)
class ChemicalPathwayRecord:
    """Observed chemical-to-pathway collective perturbation relationship from Stage 8C."""

    exposure: ChemicalExposureContext
    pathway_id: str
    pathway_name: str
    source_database: str  # "GO_BP" or "KEGG"
    database_version: str | None
    enrichment_score: float
    normalized_enrichment_score: float
    p_value: float
    adjusted_p_value: float
    pathway_size: int
    pathway_direction: str  # "POSITIVE" or "NEGATIVE"
    is_significant: bool
    leading_edge_gene_ids: tuple[str, ...]
    leading_edge_symbols: tuple[str, ...]
    all_member_gene_ids: tuple[str, ...] = ()
    all_member_symbols: tuple[str, ...] = ()
    has_complete_membership: bool = False
    ranking_metric: str = "wald_statistic"
    analysis_version: str | None = None
    scientific_limitation: str = (
        "Statistically evaluated pathway enrichment from preranked GSEA. "
        "Reflects coordinated transcriptional modulation across the gene set."
    )

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["exposure"] = self.exposure.to_dict()
        d["leading_edge_gene_ids"] = list(self.leading_edge_gene_ids)
        d["leading_edge_symbols"] = list(self.leading_edge_symbols)
        d["all_member_gene_ids"] = list(self.all_member_gene_ids)
        d["all_member_symbols"] = list(self.all_member_symbols)
        return d


@dataclass(frozen=True, slots=True)
class ChemicalGenePathwayEvidenceRecord:
    """The canonical Stage 8D evidence unit: Chemical Exposure x Observed Gene x Pathway.

    Preserves full context, original effect sizes, pathway membership, and explicit
    scientific limitations in an immutable structured record.
    """

    evidence_id: str
    study_id: str
    contrast_id: str
    independent_study_id: str | None
    compound_id: str
    compound_name: str
    treatment_condition_id: str | None
    matched_control_condition_ids: tuple[str, ...]
    dose_value: float | None
    dose_unit: str | None
    exposure_duration: float | None
    exposure_duration_unit: str | None
    developmental_stage: str | None
    model_type: str | None

    # Gene Evidence
    canonical_ensembl_gene_id: str
    gene_symbol: str | None
    base_mean: float | None
    log2_fold_change: float | None
    log2_fold_change_se: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value: float | None
    gene_response_direction: ResponseDirection
    gene_estimate_status: EstimateStatus

    # Pathway Evidence
    pathway_evidence_status: PathwayEvidenceStatus
    pathway_id: str | None
    pathway_name: str | None
    source_database: str | None
    database_version: str | None = None
    enrichment_score: float | None = None
    pathway_normalized_enrichment_score: float | None = None
    pathway_p_value: float | None = None
    pathway_adjusted_p_value: float | None = None
    pathway_size: int | None = None
    pathway_direction: str | None = None
    is_pathway_significant: bool | None = None

    # Relationship Context
    gene_membership_type: GenePathwayMembershipType = GenePathwayMembershipType.UNKNOWN_OR_UNAVAILABLE
    is_leading_edge: bool = False
    evidence_type: EvidenceType = EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION
    is_direct_target_interaction: bool = False
    causality_claim: str = "OBSERVED_TRANSCRIPTOMIC_RESPONSE"
    scientific_limitations: str = ""
    provenance_notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        ctrl_str = ";".join(self.matched_control_condition_ids) if self.matched_control_condition_ids else ""
        return {
            "evidence_id": self.evidence_id,
            "study_id": self.study_id,
            "contrast_id": self.contrast_id,
            "independent_study_id": self.independent_study_id,
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "treatment_condition_id": self.treatment_condition_id,
            "matched_control_condition_ids": ctrl_str,
            "dose_value": self.dose_value,
            "dose_unit": self.dose_unit,
            "exposure_duration": self.exposure_duration,
            "exposure_duration_unit": self.exposure_duration_unit,
            "developmental_stage": self.developmental_stage,
            "model_type": self.model_type,
            "canonical_ensembl_gene_id": self.canonical_ensembl_gene_id,
            "gene_symbol": self.gene_symbol,
            "base_mean": self.base_mean,
            "log2_fold_change": self.log2_fold_change,
            "log2_fold_change_se": self.log2_fold_change_se,
            "wald_statistic": self.wald_statistic,
            "p_value": self.p_value,
            "adjusted_p_value": self.adjusted_p_value,
            "gene_response_direction": str(self.gene_response_direction),
            "gene_estimate_status": str(self.gene_estimate_status),
            "pathway_evidence_status": str(self.pathway_evidence_status),
            "pathway_id": self.pathway_id,
            "pathway_name": self.pathway_name,
            "source_database": self.source_database,
            "database_version": self.database_version,
            "enrichment_score": self.enrichment_score,
            "pathway_normalized_enrichment_score": self.pathway_normalized_enrichment_score,
            "pathway_p_value": self.pathway_p_value,
            "pathway_adjusted_p_value": self.pathway_adjusted_p_value,
            "pathway_size": self.pathway_size,
            "pathway_direction": self.pathway_direction,
            "is_pathway_significant": self.is_pathway_significant,
            "gene_membership_type": str(self.gene_membership_type),
            "is_leading_edge": self.is_leading_edge,
            "evidence_type": str(self.evidence_type),
            "is_direct_target_interaction": self.is_direct_target_interaction,
            "causality_claim": self.causality_claim,
            "scientific_limitations": self.scientific_limitations,
            "provenance_notes": self.provenance_notes,
        }


@dataclass(frozen=True, slots=True)
class ChemicalGenePathwayEvidenceDataset:
    """Governed aggregate container for Stage 8D evidence relationships."""

    dataset_id: str
    records: tuple[ChemicalGenePathwayEvidenceRecord, ...]
    chemical_gene_records: tuple[ChemicalGeneRecord, ...]
    chemical_pathway_records: tuple[ChemicalPathwayRecord, ...]
    created_at: str
    analysis_version: str = "1.0.0"

    @property
    def total_evidence_records(self) -> int:
        return len(self.records)

    @property
    def total_chemical_gene_records(self) -> int:
        return len(self.chemical_gene_records)

    @property
    def total_chemical_pathway_records(self) -> int:
        return len(self.chemical_pathway_records)

    @property
    def unique_compounds(self) -> set[str]:
        return {r.compound_id for r in self.records if r.compound_id}

    @property
    def unique_contrasts(self) -> set[str]:
        return {r.contrast_id for r in self.records if r.contrast_id}

    @property
    def unique_genes(self) -> set[str]:
        return {r.canonical_ensembl_gene_id for r in self.records if r.canonical_ensembl_gene_id}

    @property
    def unique_pathways(self) -> set[str]:
        return {r.pathway_id for r in self.records if r.pathway_id}

    def get_records_for_gene(self, gene_id: str) -> tuple[ChemicalGenePathwayEvidenceRecord, ...]:
        return tuple(r for r in self.records if r.canonical_ensembl_gene_id == gene_id)

    def get_records_for_pathway(self, pathway_id: str) -> tuple[ChemicalGenePathwayEvidenceRecord, ...]:
        return tuple(r for r in self.records if r.pathway_id == pathway_id)

    def get_records_for_contrast(self, contrast_id: str) -> tuple[ChemicalGenePathwayEvidenceRecord, ...]:
        return tuple(r for r in self.records if r.contrast_id == contrast_id)

    def to_dict(self) -> dict[str, Any]:
        return {
            "dataset_id": self.dataset_id,
            "created_at": self.created_at,
            "analysis_version": self.analysis_version,
            "total_evidence_records": self.total_evidence_records,
            "total_chemical_gene_records": self.total_chemical_gene_records,
            "total_chemical_pathway_records": self.total_chemical_pathway_records,
            "unique_compounds_count": len(self.unique_compounds),
            "unique_contrasts_count": len(self.unique_contrasts),
            "unique_genes_count": len(self.unique_genes),
            "unique_pathways_count": len(self.unique_pathways),
            "records": [r.to_dict() for r in self.records],
        }
