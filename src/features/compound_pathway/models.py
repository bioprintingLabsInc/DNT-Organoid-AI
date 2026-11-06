"""Stage 8D: Chemical–Gene–Pathway Evidence Relationships — Data Models.

Governed Scientific Invariants:
1. Tripartite Traceability: Connects:
       Chemical Exposure Context -> Observed Gene Response -> Associated Pathway Membership
   while preserving the exact upstream statistics and experimental parameters.
2. No Causality Assumption: Distinguishes observed transcriptomic perturbation from
   experimentally proven direct chemical-target binding or receptor interaction.
3. No Direction Transference: A gene's individual log2FC sign does NOT infer overall pathway
   activation or inhibition without formal pathway-level scoring.
4. Universal Modality Independence: Agnostic to chemical, dose, duration, organoid age,
   or RNA-seq platform. Zero hardcoded validation dataset constants.
5. Strict Leakage Prevention: Contains zero Stage 9 DNT outcome labels or supervisory signals.
"""

from dataclasses import asdict, dataclass
from enum import StrEnum
import math
from typing import Any

from src.features.cross_study.models import (
    EstimateStatus,
    ResponseDirection,
)


class EvidenceType(StrEnum):
    """Scientific nature of the observed evidence connection."""

    OBSERVED_TRANSCRIPTOMIC_PERTURBATION = "OBSERVED_TRANSCRIPTOMIC_PERTURBATION"
    PATHWAY_ENRICHMENT = "PATHWAY_ENRICHMENT"
    DIRECT_TARGET_INTERACTION = "DIRECT_TARGET_INTERACTION"
    PREDICTED_ASSOCIATION = "PREDICTED_ASSOCIATION"


class PathwayEvidenceStatus(StrEnum):
    """Availability and qualification status of pathway enrichment data for a contrast."""

    AVAILABLE = "AVAILABLE"
    NOT_AVAILABLE = "NOT_AVAILABLE"
    EVALUATION_INELIGIBLE = "EVALUATION_INELIGIBLE"
    INSUFFICIENT_GENES = "INSUFFICIENT_GENES"


class GenePathwayMembershipType(StrEnum):
    """Nature of gene participation in a biological pathway."""

    LEADING_EDGE_MEMBER = "LEADING_EDGE_MEMBER"
    PATHWAY_ANNOTATED_MEMBER = "PATHWAY_ANNOTATED_MEMBER"
    UNANNOTATED = "UNANNOTATED"


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
    biological_source_ids: tuple[str, ...]
    organoid_context_ids: tuple[str, ...]
    study_id: str
    study_family_id: str
    independent_study_id: str
    experiment_id: str
    contrast_id: str
    treatment_condition_id: str
    matched_control_condition_ids: tuple[str, ...]
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]
    treatment_replicate_count: int
    control_replicate_count: int

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
    base_mean: float
    log2_fold_change: float | None
    log2_fold_change_se: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value: float | None
    response_direction: ResponseDirection
    estimate_status: EstimateStatus
    evidence_type: EvidenceType = EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION
    is_direct_target_interaction: bool = False
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
    ranking_metric: str
    analysis_version: str
    scientific_limitation: str = (
        "Statistically evaluated pathway enrichment from preranked GSEA. "
        "Reflects coordinated transcriptional modulation across the gene set."
    )

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["exposure"] = self.exposure.to_dict()
        d["leading_edge_gene_ids"] = list(self.leading_edge_gene_ids)
        d["leading_edge_symbols"] = list(self.leading_edge_symbols)
        return d


@dataclass(frozen=True, slots=True)
class ChemicalGenePathwayEvidenceRecord:
    """The canonical Stage 8D evidence unit: Chemical Exposure x Observed Gene x Pathway.

    Preserves full context, original effect sizes, pathway membership, and explicit
    scientific limitations in an immutable structured record.
    """

    evidence_id: str
    contrast_id: str
    study_id: str
    independent_study_id: str
    compound_id: str
    compound_name: str
    dose_value: float | None
    dose_unit: str | None
    exposure_duration: float | None
    exposure_duration_unit: str | None
    developmental_stage: str | None
    model_type: str | None

    # Gene Evidence
    canonical_ensembl_gene_id: str
    gene_symbol: str | None
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
    database_version: str | None
    pathway_normalized_enrichment_score: float | None
    pathway_p_value: float | None
    pathway_adjusted_p_value: float | None
    pathway_direction: str | None
    is_pathway_significant: bool | None

    # Relationship Context
    gene_membership_type: GenePathwayMembershipType
    is_leading_edge: bool
    evidence_type: EvidenceType
    is_direct_target_interaction: bool
    scientific_limitations: str
    provenance_notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "evidence_id": self.evidence_id,
            "contrast_id": self.contrast_id,
            "study_id": self.study_id,
            "independent_study_id": self.independent_study_id,
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "dose_value": self.dose_value,
            "dose_unit": self.dose_unit,
            "exposure_duration": self.exposure_duration,
            "exposure_duration_unit": self.exposure_duration_unit,
            "developmental_stage": self.developmental_stage,
            "model_type": self.model_type,
            "canonical_ensembl_gene_id": self.canonical_ensembl_gene_id,
            "gene_symbol": self.gene_symbol,
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
            "pathway_normalized_enrichment_score": self.pathway_normalized_enrichment_score,
            "pathway_p_value": self.pathway_p_value,
            "pathway_adjusted_p_value": self.pathway_adjusted_p_value,
            "pathway_direction": self.pathway_direction,
            "is_pathway_significant": self.is_pathway_significant,
            "gene_membership_type": str(self.gene_membership_type),
            "is_leading_edge": self.is_leading_edge,
            "evidence_type": str(self.evidence_type),
            "is_direct_target_interaction": self.is_direct_target_interaction,
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
