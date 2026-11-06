"""Tests for Stage 8B Network Topology Features and Scientific Interpretation Guardrails.

Verifies:
1. Computation of graph topology metrics (degree, clustering coefficient, neighborhood size, percentile rank).
2. Enforcement of REFERENCE_DEGREE_BIAS_DISCLAIMER on all records.
3. Graceful handling of unannotated or isolated genes.
4. Chemical-responsive induced subnetwork extraction.
5. Community / functional module detection via greedy modularity.
6. Guardrails preventing unsupported causal, binding, or biomarker claims.
"""

import pytest
from src.features.coresponse.models import (
    BiologicalNetwork,
    EvidenceCategory,
    GeneNetworkEdge,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
)
from src.features.coresponse.topology import (
    REFERENCE_DEGREE_BIAS_DISCLAIMER,
    build_networkx_graph,
    compute_gene_network_topology,
    detect_network_modules,
    extract_chemical_responsive_subnetwork,
)
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


def _build_synthetic_triangle_network() -> BiologicalNetwork:
    """Build a network with 3 interconnected genes (triangle: G1-G2, G2-G3, G3-G1) plus a branch G3-G4."""
    edges = (
        _edge(
            gene_a_canonical_id="ENSG000001",
            gene_b_canonical_id="ENSG000002",
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            confidence_score=0.9,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ),
        _edge(
            gene_a_canonical_id="ENSG000002",
            gene_b_canonical_id="ENSG000003",
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            confidence_score=0.85,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ),
        _edge(
            gene_a_canonical_id="ENSG000001",
            gene_b_canonical_id="ENSG000003",
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            confidence_score=0.92,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ),
        _edge(
            gene_a_canonical_id="ENSG000003",
            gene_b_canonical_id="ENSG000004",
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.FUNCTIONAL_ASSOCIATION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            confidence_score=0.7,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ),
    )
    return BiologicalNetwork(
        network_id="synth_net_triangle",
        network_name="Synthetic Network for Topology Tests",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        edges=edges,
        source_database="Synthetic",
        database_version="1.0",
    )


def test_topology_metrics_calculation_and_degree_bias_disclaimer() -> None:
    """Verify exact degree, clustering coefficient, percentile rank, and mandatory bias note."""
    net = _build_synthetic_triangle_network()
    symbols = {
        "ENSG000001": "GENE1",
        "ENSG000002": "GENE2",
        "ENSG000003": "GENE3",
        "ENSG000004": "GENE4",
    }

    metrics = compute_gene_network_topology(
        network=net,
        target_genes=["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004", "ENSG_ABSENT"],
        gene_symbol_map=symbols,
    )

    # G3 is connected to G1, G2, G4 -> degree 3
    assert metrics["ENSG000003"].degree == 3
    assert metrics["ENSG000003"].neighborhood_size == 3
    assert not metrics["ENSG000003"].is_isolated
    assert metrics["ENSG000003"].gene_symbol == "GENE3"

    # G1 is connected to G2 and G3 (which are connected) -> clustering coeff = 1.0
    assert metrics["ENSG000001"].degree == 2
    assert pytest.approx(metrics["ENSG000001"].clustering_coefficient, 0.01) == 1.0

    # G4 is connected only to G3 -> degree 1, clustering = 0.0
    assert metrics["ENSG000004"].degree == 1
    assert metrics["ENSG000004"].clustering_coefficient == 0.0

    # Unconnected/absent gene
    assert metrics["ENSG_ABSENT"].degree == 0
    assert metrics["ENSG_ABSENT"].is_isolated is True
    assert metrics["ENSG_ABSENT"].clustering_coefficient == 0.0
    assert metrics["ENSG_ABSENT"].degree_percentile == 0.0

    # All evaluated genes must include the mandatory bias disclaimer
    for g_id, m in metrics.items():
        assert m.reference_degree_bias_note == REFERENCE_DEGREE_BIAS_DISCLAIMER
        assert "causal" in m.reference_degree_bias_note.lower()
        assert "biomarker" in m.reference_degree_bias_note.lower()


def test_extract_chemical_responsive_subnetwork() -> None:
    """Verify subnetwork extraction isolates edges connecting responsive genes."""
    net = _build_synthetic_triangle_network()
    # Assume only G1 and G2 were responsive to the chemical
    sub_net = extract_chemical_responsive_subnetwork(
        network=net,
        responsive_genes=["ENSG000001", "ENSG000002"],
        subnetwork_id="chem_resp_sub",
    )

    assert sub_net.network_id == "chem_resp_sub"
    assert len(sub_net.edges) == 1
    edge = sub_net.edges[0]
    assert {edge.gene_a_canonical_id, edge.gene_b_canonical_id} == {"ENSG000001", "ENSG000002"}
    assert "Induced subnetwork connecting 2 chemical-responsive genes" in sub_net.scientific_limitations
    assert "NOT imply direct chemical target interaction" in sub_net.scientific_limitations


def test_detect_network_modules_via_modularity() -> None:
    """Verify community detection partitions nodes into discrete modules."""
    # Build two distinct cliques: {G1, G2, G3} and {G10, G11, G12}
    edges = []
    # Module A: G1, G2, G3 fully connected
    for g_a, g_b in [("G1", "G2"), ("G2", "G3"), ("G1", "G3")]:
        edges.append(_edge(
            gene_a_canonical_id=g_a,
            gene_b_canonical_id=g_b,
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ))
    # Module B: G10, G11, G12 fully connected
    for g_a, g_b in [("G10", "G11"), ("G11", "G12"), ("G10", "G12")]:
        edges.append(_edge(
            gene_a_canonical_id=g_a,
            gene_b_canonical_id=g_b,
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
            evidence_source="STRING",
            evidence_version="v12.0",
            evidence_category=EvidenceCategory.DATABASE_CURATED,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        ))

    net = BiologicalNetwork(
        network_id="two_modules",
        network_name="Two Modules Network",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        edges=tuple(edges),
        source_database="Synthetic",
    )

    modules = detect_network_modules(net, min_module_size=3)
    assert len(modules) == 2
    mod_members = [set(m) for m in modules.values()]
    assert {"G1", "G2", "G3"} in mod_members
    assert {"G10", "G11", "G12"} in mod_members
