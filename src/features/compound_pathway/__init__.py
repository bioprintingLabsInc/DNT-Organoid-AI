"""Stage 8D: Chemical–Gene–Pathway Evidence Relationships Package.

Synthesizes structured, traceable relationships connecting:
    Chemical Exposure -> Observed Gene Perturbation -> Pathway Context
preserving all upstream statistics and maintaining universal chemical/dataset independence.
"""

from src.features.compound_pathway.builder import (
    build_chemical_exposure_context,
    build_chemical_gene_records,
    build_chemical_pathway_records,
    link_chemical_gene_pathway_evidence,
)
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
from src.features.compound_pathway.qc import (
    EvidenceQCSummary,
    audit_evidence_dataset_qc,
)
from src.features.compound_pathway.runner import (
    build_stage8d_evidence_pipeline,
    load_pathway_results_from_csv,
    load_registry_from_csv,
)
from src.features.compound_pathway.serialization import (
    CANONICAL_EVIDENCE_COLUMNS,
    export_evidence_to_csv,
    export_evidence_to_dataframe,
)

__all__ = [
    "CANONICAL_EVIDENCE_COLUMNS",
    "ChemicalExposureContext",
    "ChemicalGenePathwayEvidenceDataset",
    "ChemicalGenePathwayEvidenceRecord",
    "ChemicalGeneRecord",
    "ChemicalPathwayRecord",
    "EvidenceQCSummary",
    "EvidenceType",
    "GenePathwayMembershipType",
    "PathwayEvidenceStatus",
    "audit_evidence_dataset_qc",
    "build_chemical_exposure_context",
    "build_chemical_gene_records",
    "build_chemical_pathway_records",
    "build_stage8d_evidence_pipeline",
    "export_evidence_to_csv",
    "export_evidence_to_dataframe",
    "link_chemical_gene_pathway_evidence",
    "load_pathway_results_from_csv",
    "load_registry_from_csv",
]
