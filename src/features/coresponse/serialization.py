"""Stage 8B: Serialization and Governed Export for Network Features.

Provides deterministic conversion and export of:
1. Canonical Biological Network Edges to DataFrame and CSV.
2. Gene Network Topology Metrics to DataFrame and CSV.
"""

from pathlib import Path
from typing import Sequence
import pandas as pd

from src.features.coresponse.models import (
    BiologicalNetwork,
    GeneNetworkTopologyMetrics,
)


CANONICAL_NETWORK_EDGE_COLUMNS = [
    "edge_id",
    "gene_a_canonical_id",
    "gene_b_canonical_id",
    "gene_a_symbol",
    "gene_b_symbol",
    "relationship_type",
    "interaction_nature",
    "evidence_source",
    "evidence_version",
    "evidence_category",
    "confidence_score",
    "statistical_metric_name",
    "statistical_estimate",
    "p_value",
    "adjusted_p_value",
    "sample_size",
    "condition_count",
    "study_context",
    "compound_context",
    "partition_assignment",
    "qc_status",
    "is_direct_binding",
    "is_chemical_specific",
    "causality_claim",
    "scientific_limitations",
    "provenance_notes",
]

CANONICAL_TOPOLOGY_METRICS_COLUMNS = [
    "canonical_ensembl_gene_id",
    "gene_symbol",
    "network_id",
    "degree",
    "clustering_coefficient",
    "neighborhood_size",
    "is_isolated",
    "degree_percentile",
    "reference_degree_bias_note",
]


def export_network_to_dataframe(network: BiologicalNetwork) -> pd.DataFrame:
    """Convert BiologicalNetwork edges to a canonical DataFrame."""
    rows = [e.to_dict() for e in network.edges]
    df = pd.DataFrame(rows, columns=CANONICAL_NETWORK_EDGE_COLUMNS)
    if not df.empty:
        df = df.sort_values(
            by=["gene_a_canonical_id", "gene_b_canonical_id", "relationship_type"],
            ascending=True,
            na_position="last",
        )
    return df


def export_network_to_csv(network: BiologicalNetwork, output_path: str | Path) -> Path:
    """Export BiologicalNetwork edges deterministically to CSV."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = export_network_to_dataframe(network)
    df.to_csv(path, index=False, float_format="%.6f")
    return path


def export_topology_metrics_to_dataframe(
    metrics_map: dict[str, GeneNetworkTopologyMetrics],
) -> pd.DataFrame:
    """Convert gene topology metrics dictionary to a canonical DataFrame."""
    rows = [m.to_dict() for m in metrics_map.values()]
    df = pd.DataFrame(rows, columns=CANONICAL_TOPOLOGY_METRICS_COLUMNS)
    if not df.empty:
        df = df.sort_values(by=["canonical_ensembl_gene_id"], ascending=True)
    return df


def export_topology_metrics_to_csv(
    metrics_map: dict[str, GeneNetworkTopologyMetrics],
    output_path: str | Path,
) -> Path:
    """Export gene topology metrics deterministically to CSV."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = export_topology_metrics_to_dataframe(metrics_map)
    df.to_csv(path, index=False, float_format="%.6f")
    return path
