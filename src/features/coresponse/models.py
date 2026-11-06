"""Stage 8B: Data Models and Schemas for Gene Co-Response and Biological Network Features.

Defines immutable, typed dataclasses, controlled vocabularies, and governed schemas
strictly distinguishing:
1. Established biological reference networks (independent of chemical exposure).
2. Observed within-study gene co-expression (sample-level continuous expression).
3. Cross-condition gene co-response (coordinated response across qualified perturbation contrasts).
4. Data-adaptive prediction features (strictly partitioned to prevent AI/ML leakage).
"""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import StrEnum
import hashlib
from typing import Any


class NetworkRelationshipType(StrEnum):
    """Governed taxonomy of biological network evidence types."""

    ESTABLISHED_REFERENCE_INTERACTION = "ESTABLISHED_REFERENCE_INTERACTION"
    WITHIN_STUDY_CO_EXPRESSION = "WITHIN_STUDY_CO_EXPRESSION"
    CROSS_CONDITION_CO_RESPONSE = "CROSS_CONDITION_CO_RESPONSE"
    DATA_ADAPTIVE_PREDICTION_FEATURE = "DATA_ADAPTIVE_PREDICTION_FEATURE"


class InteractionNature(StrEnum):
    """Specific biological or empirical nature of the gene-gene relationship."""

    PHYSICAL_INTERACTION = "PHYSICAL_INTERACTION"
    FUNCTIONAL_ASSOCIATION = "FUNCTIONAL_ASSOCIATION"
    PREDICTED_INTERACTION = "PREDICTED_INTERACTION"
    STATISTICAL_CO_EXPRESSION = "STATISTICAL_CO_EXPRESSION"
    STATISTICAL_CO_RESPONSE = "STATISTICAL_CO_RESPONSE"
    REGULATORY_INTERACTION = "REGULATORY_INTERACTION"
    DESCRIPTIVE_CONCORDANCE = "DESCRIPTIVE_CONCORDANCE"
    UNKNOWN_OR_UNSPECIFIED = "UNKNOWN_OR_UNSPECIFIED"


class EvidenceCategory(StrEnum):
    """Source authority or derivation mechanism of relationship evidence."""

    EXPERIMENTAL_ASSAY = "EXPERIMENTAL_ASSAY"
    DATABASE_CURATED = "DATABASE_CURATED"
    TEXT_MINING = "TEXT_MINING"
    COMPUTATIONAL_PREDICTION = "COMPUTATIONAL_PREDICTION"
    EMPIRICAL_COEXPRESSION = "EMPIRICAL_COEXPRESSION"
    EMPIRICAL_CORESPONSE = "EMPIRICAL_CORESPONSE"
    DESCRIPTIVE_CONCORDANCE = "DESCRIPTIVE_CONCORDANCE"
    COMBINED_SCORE = "COMBINED_SCORE"


class PartitionAssignment(StrEnum):
    """AI/ML partitioning safety governance for network evidence."""

    UNRESTRICTED_REFERENCE = "UNRESTRICTED_REFERENCE"
    DESCRIPTIVE_EXPLORATION = "DESCRIPTIVE_EXPLORATION"
    TRAINING_ONLY = "TRAINING_ONLY"
    HELD_OUT_EVALUATION = "HELD_OUT_EVALUATION"


class CoresponseQCStatus(StrEnum):
    """Quality control and statistical eligibility status for network inference."""

    VALID = "VALID"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    INSUFFICIENT_SAMPLES = "INSUFFICIENT_SAMPLES"
    INSUFFICIENT_CONDITIONS = "INSUFFICIENT_CONDITIONS"
    INSUFFICIENT_CHEMICALS = "INSUFFICIENT_CHEMICALS"
    UNANNOTATED = "UNANNOTATED"
    LEAKAGE_DETECTED = "LEAKAGE_DETECTED"


def derive_edge_id(gene_a: str, gene_b: str, relationship_type: str, evidence_source: str) -> str:
    """Derive deterministic canonical edge identifier sorted alphabetically by gene ID."""
    g1, g2 = sorted([gene_a.strip().upper(), gene_b.strip().upper()])
    raw_key = f"{g1}::{g2}::{relationship_type.strip()}::{evidence_source.strip()}"
    digest = hashlib.sha256(raw_key.encode("utf-8")).hexdigest()[:16]
    return f"edge_{digest}"


@dataclass(frozen=True, slots=True)
class GeneNetworkEdge:
    """Canonical edge record connecting two genes in a biological or empirical network."""

    edge_id: str
    gene_a_canonical_id: str
    gene_b_canonical_id: str
    gene_a_symbol: str | None = None
    gene_b_symbol: str | None = None
    relationship_type: NetworkRelationshipType = NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION
    interaction_nature: InteractionNature = InteractionNature.UNKNOWN_OR_UNSPECIFIED
    evidence_source: str = ""
    evidence_version: str | None = None
    evidence_category: EvidenceCategory = EvidenceCategory.DATABASE_CURATED
    confidence_score: float | None = None
    statistical_metric_name: str | None = None
    statistical_estimate: float | None = None
    p_value: float | None = None
    adjusted_p_value: float | None = None
    sample_size: int | None = None
    condition_count: int | None = None
    study_context: tuple[str, ...] = ()
    compound_context: tuple[str, ...] = ()
    partition_assignment: PartitionAssignment = PartitionAssignment.UNRESTRICTED_REFERENCE
    qc_status: CoresponseQCStatus = CoresponseQCStatus.VALID
    is_direct_binding: bool = False
    is_chemical_specific: bool = False
    causality_claim: str = "NON_CAUSAL_ASSOCIATION"
    scientific_limitations: str = ""
    provenance_notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Convert edge to dictionary for governed serialization."""
        return {
            "edge_id": self.edge_id,
            "gene_a_canonical_id": self.gene_a_canonical_id,
            "gene_b_canonical_id": self.gene_b_canonical_id,
            "gene_a_symbol": self.gene_a_symbol,
            "gene_b_symbol": self.gene_b_symbol,
            "relationship_type": str(self.relationship_type),
            "interaction_nature": str(self.interaction_nature),
            "evidence_source": self.evidence_source,
            "evidence_version": self.evidence_version,
            "evidence_category": str(self.evidence_category),
            "confidence_score": self.confidence_score,
            "statistical_metric_name": self.statistical_metric_name,
            "statistical_estimate": self.statistical_estimate,
            "p_value": self.p_value,
            "adjusted_p_value": self.adjusted_p_value,
            "sample_size": self.sample_size,
            "condition_count": self.condition_count,
            "study_context": ";".join(self.study_context) if self.study_context else None,
            "compound_context": ";".join(self.compound_context) if self.compound_context else None,
            "partition_assignment": str(self.partition_assignment),
            "qc_status": str(self.qc_status),
            "is_direct_binding": self.is_direct_binding,
            "is_chemical_specific": self.is_chemical_specific,
            "causality_claim": self.causality_claim,
            "scientific_limitations": self.scientific_limitations,
            "provenance_notes": self.provenance_notes,
        }


@dataclass(frozen=True, slots=True)
class BiologicalNetwork:
    """An immutable, versioned network collection of gene-gene edges."""

    network_id: str
    network_name: str
    relationship_type: NetworkRelationshipType
    edges: tuple[GeneNetworkEdge, ...]
    species: str = "Homo sapiens"
    taxon_id: int = 9606
    source_database: str | None = None
    database_version: str | None = None
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    qc_status: CoresponseQCStatus = CoresponseQCStatus.VALID
    scientific_limitations: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def edge_count(self) -> int:
        return len(self.edges)

    @property
    def unique_genes(self) -> set[str]:
        genes: set[str] = set()
        for e in self.edges:
            genes.add(e.gene_a_canonical_id)
            genes.add(e.gene_b_canonical_id)
        return genes

    def get_edges_for_gene(self, gene_id: str) -> tuple[GeneNetworkEdge, ...]:
        gid = gene_id.strip().upper()
        return tuple(
            e for e in self.edges
            if e.gene_a_canonical_id.upper() == gid or e.gene_b_canonical_id.upper() == gid
        )


@dataclass(frozen=True, slots=True)
class GeneNetworkTopologyMetrics:
    """Graph topology and connectivity metrics for a single gene."""

    canonical_ensembl_gene_id: str
    gene_symbol: str | None
    network_id: str
    degree: int
    clustering_coefficient: float
    neighborhood_size: int
    is_isolated: bool
    degree_percentile: float
    reference_degree_bias_note: str = (
        "Gene degree reflects reference literature investigation density as well as "
        "biological connectivity. High degree does NOT establish causal criticality or biomarker status."
    )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CoresponseEvaluationSummary:
    """Summary of co-response or co-expression evaluation and statistical gating."""

    status: CoresponseQCStatus
    is_eligible: bool
    evaluated_genes_count: int
    evaluated_conditions_count: int
    evaluated_chemicals_count: int
    evaluated_samples_count: int
    inferred_edges_count: int
    reasons: tuple[str, ...] = ()
    scientific_limitations: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "is_eligible": self.is_eligible,
            "evaluated_genes_count": self.evaluated_genes_count,
            "evaluated_conditions_count": self.evaluated_conditions_count,
            "evaluated_chemicals_count": self.evaluated_chemicals_count,
            "evaluated_samples_count": self.evaluated_samples_count,
            "inferred_edges_count": self.inferred_edges_count,
            "reasons": list(self.reasons),
            "scientific_limitations": self.scientific_limitations,
        }
