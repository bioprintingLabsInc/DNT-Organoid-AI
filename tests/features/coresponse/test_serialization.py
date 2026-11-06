"""Tests for Stage 8B Serialization and Governed Output Export.

Verifies:
1. Exact compliance with the 26-column canonical network edge schema.
2. Exact compliance with the 9-column canonical topology metrics schema.
3. Deterministic alphabetical sorting by canonical Ensembl gene IDs.
4. Correct float formatting and round-trip CSV serialization.
"""

from pathlib import Path
import pandas as pd
import pytest

from src.features.coresponse.models import (
    BiologicalNetwork,
    EvidenceCategory,
    GeneNetworkEdge,
    GeneNetworkTopologyMetrics,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
)
from src.features.coresponse.serialization import (
    CANONICAL_NETWORK_EDGE_COLUMNS,
    CANONICAL_TOPOLOGY_METRICS_COLUMNS,
    export_network_to_csv,
    export_network_to_dataframe,
    export_topology_metrics_to_csv,
    export_topology_metrics_to_dataframe,
)
from src.features.coresponse.topology import REFERENCE_DEGREE_BIAS_DISCLAIMER
from src.features.coresponse.models import derive_edge_id


def _edge(**kwargs) -> GeneNetworkEdge:
    """Construct an edge with its deterministic canonical edge_id derived by production code."""
    kwargs["edge_id"] = derive_edge_id(
        kwargs["gene_a_canonical_id"],
        kwargs["gene_b_canonical_id"],
        str(kwargs["relationship_type"]),
        kwargs["evidence_source"],
    )
    return GeneNetworkEdge(**kwargs)


def test_export_network_to_dataframe_schema_and_sorting() -> None:
    """Verify DataFrame export adheres strictly to the 26-column schema and sorts deterministically."""
    e1 = _edge(
        gene_a_canonical_id="ENSG000002",
        gene_b_canonical_id="ENSG000003",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
        evidence_source="STRING",
        evidence_version="v12.0",
        evidence_category=EvidenceCategory.DATABASE_CURATED,
        confidence_score=0.9,
    )
    e2 = _edge(
        gene_a_canonical_id="ENSG000001",
        gene_b_canonical_id="ENSG000004",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        interaction_nature=InteractionNature.FUNCTIONAL_ASSOCIATION,
        evidence_source="STRING",
        evidence_version="v12.0",
        evidence_category=EvidenceCategory.DATABASE_CURATED,
        confidence_score=0.8,
    )

    net = BiologicalNetwork(
        network_id="test_export_net",
        network_name="Export Network",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        edges=(e1, e2),
        source_database="STRING",
    )

    df = export_network_to_dataframe(net)

    # Verify columns exactly match the 26-column schema
    assert list(df.columns) == CANONICAL_NETWORK_EDGE_COLUMNS
    assert len(df) == 2

    # Deterministic sorting: ENSG000001 must appear before ENSG000002
    assert df.iloc[0]["gene_a_canonical_id"] == "ENSG000001"
    assert df.iloc[1]["gene_a_canonical_id"] == "ENSG000002"


def test_export_network_to_csv(tmp_path: Path) -> None:
    """Verify CSV export writes valid, re-readable CSV with correct formatting."""
    e = _edge(
        gene_a_canonical_id="ENSG000001",
        gene_b_canonical_id="ENSG000002",
        relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
        interaction_nature=InteractionNature.STATISTICAL_CO_EXPRESSION,
        evidence_source="Stage8B_Coexpression",
        evidence_version="1.0.0",
        evidence_category=EvidenceCategory.EMPIRICAL_COEXPRESSION,
        statistical_metric_name="spearman_rho",
        statistical_estimate=0.876543,
        p_value=0.0001234,
        adjusted_p_value=0.0004567,
        sample_size=16,
        study_context=("GSE999999",),
        partition_assignment=PartitionAssignment.DESCRIPTIVE_EXPLORATION,
    )
    net = BiologicalNetwork(
        network_id="csv_net",
        network_name="CSV Net",
        relationship_type=NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
        edges=(e,),
        source_database="Empirical",
    )

    csv_file = tmp_path / "network_edges.csv"
    res_path = export_network_to_csv(net, csv_file)

    assert res_path.exists()
    re_df = pd.read_csv(res_path)
    assert list(re_df.columns) == CANONICAL_NETWORK_EDGE_COLUMNS
    assert len(re_df) == 1
    assert re_df.iloc[0]["sample_size"] == 16
    assert pytest.approx(re_df.iloc[0]["statistical_estimate"], 1e-4) == 0.876543


def test_export_topology_metrics_to_dataframe_and_csv(tmp_path: Path) -> None:
    """Verify topology metrics DataFrame and CSV export strictly matches 9-column schema."""
    m1 = GeneNetworkTopologyMetrics(
        canonical_ensembl_gene_id="ENSG000002",
        gene_symbol="GENE2",
        network_id="net_1",
        degree=5,
        clustering_coefficient=0.25,
        neighborhood_size=5,
        is_isolated=False,
        degree_percentile=75.0,
        reference_degree_bias_note=REFERENCE_DEGREE_BIAS_DISCLAIMER,
    )
    m2 = GeneNetworkTopologyMetrics(
        canonical_ensembl_gene_id="ENSG000001",
        gene_symbol="GENE1",
        network_id="net_1",
        degree=10,
        clustering_coefficient=0.5,
        neighborhood_size=10,
        is_isolated=False,
        degree_percentile=95.0,
        reference_degree_bias_note=REFERENCE_DEGREE_BIAS_DISCLAIMER,
    )

    metrics_map = {"ENSG000002": m1, "ENSG000001": m2}
    df = export_topology_metrics_to_dataframe(metrics_map)

    # Schema verification
    assert list(df.columns) == CANONICAL_TOPOLOGY_METRICS_COLUMNS
    assert len(df) == 2
    # Sorting verification: ENSG000001 before ENSG000002
    assert df.iloc[0]["canonical_ensembl_gene_id"] == "ENSG000001"

    # CSV export
    csv_file = tmp_path / "topology_metrics.csv"
    export_topology_metrics_to_csv(metrics_map, csv_file)
    assert csv_file.exists()

    re_df = pd.read_csv(csv_file)
    assert list(re_df.columns) == CANONICAL_TOPOLOGY_METRICS_COLUMNS
    assert re_df.iloc[0]["degree"] == 10
    assert re_df.iloc[0]["reference_degree_bias_note"] == REFERENCE_DEGREE_BIAS_DISCLAIMER
