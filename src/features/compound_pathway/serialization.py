"""Stage 8D: Serialization and Export for Chemical–Gene–Pathway Evidence Relationships."""

from pathlib import Path
from typing import Any
import pandas as pd

from src.features.compound_pathway.models import ChemicalGenePathwayEvidenceDataset


CANONICAL_EVIDENCE_COLUMNS = [
    "evidence_id",
    "contrast_id",
    "study_id",
    "independent_study_id",
    "compound_id",
    "compound_name",
    "dose_value",
    "dose_unit",
    "exposure_duration",
    "exposure_duration_unit",
    "developmental_stage",
    "model_type",
    "canonical_ensembl_gene_id",
    "gene_symbol",
    "log2_fold_change",
    "log2_fold_change_se",
    "wald_statistic",
    "p_value",
    "adjusted_p_value",
    "gene_response_direction",
    "gene_estimate_status",
    "pathway_evidence_status",
    "pathway_id",
    "pathway_name",
    "source_database",
    "database_version",
    "pathway_normalized_enrichment_score",
    "pathway_p_value",
    "pathway_adjusted_p_value",
    "pathway_direction",
    "is_pathway_significant",
    "gene_membership_type",
    "is_leading_edge",
    "evidence_type",
    "is_direct_target_interaction",
    "scientific_limitations",
    "provenance_notes",
]


def export_evidence_to_dataframe(dataset: ChemicalGenePathwayEvidenceDataset) -> pd.DataFrame:
    """Convert ChemicalGenePathwayEvidenceDataset to a canonical DataFrame."""
    rows = [r.to_dict() for r in dataset.records]
    df = pd.DataFrame(rows, columns=CANONICAL_EVIDENCE_COLUMNS)
    if not df.empty:
        # Deterministic sorting by contrast, gene, and pathway
        df = df.sort_values(
            by=["contrast_id", "canonical_ensembl_gene_id", "pathway_id"],
            ascending=True,
            na_position="last",
        ).reset_index(drop=True)
    return df


def export_evidence_to_csv(
    dataset: ChemicalGenePathwayEvidenceDataset,
    output_path: str | Path,
) -> Path:
    """Export evidence dataset to a deterministic CSV file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = export_evidence_to_dataframe(dataset)
    df.to_csv(path, index=False)
    return path
