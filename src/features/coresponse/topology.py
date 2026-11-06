"""Stage 8B: Network Topology Feature Extraction and Scientific Interpretation.

Computes reproducible graph topology metrics using NetworkX while enforcing:
1. Reference-network degree bias transparency (hub genes like TP53 or EGFR have elevated degrees
   partly due to literature investigation density, not solely biological primacy).
2. Explicit guardrails preventing unsupported claims of causal regulatory control, direct target
   binding, or biomarker status.
3. Subnetwork extraction and community/module identification among chemical-perturbed genes.
"""

from collections.abc import Iterable, Sequence
from typing import Any
import networkx as nx
import numpy as np

from src.features.coresponse.models import (
    BiologicalNetwork,
    GeneNetworkTopologyMetrics,
)


REFERENCE_DEGREE_BIAS_DISCLAIMER = (
    "Gene degree reflects reference literature investigation density as well as "
    "biological connectivity. High degree does NOT establish causal criticality, "
    "direct chemical binding, or validated DNT biomarker status."
)


def build_networkx_graph(network: BiologicalNetwork) -> nx.Graph:
    """Construct an undirected NetworkX Graph from a BiologicalNetwork."""
    G = nx.Graph()
    for edge in network.edges:
        g1 = edge.gene_a_canonical_id
        g2 = edge.gene_b_canonical_id
        weight = edge.confidence_score or edge.statistical_estimate or 1.0
        G.add_edge(
            g1,
            g2,
            weight=abs(float(weight)),
            edge_id=edge.edge_id,
            relationship_type=str(edge.relationship_type),
            interaction_nature=str(edge.interaction_nature),
        )
    return G


def compute_gene_network_topology(
    network: BiologicalNetwork,
    target_genes: Iterable[str] | None = None,
    gene_symbol_map: dict[str, str] | None = None,
) -> dict[str, GeneNetworkTopologyMetrics]:
    """Compute network topology metrics for genes within a biological network.

    Args:
        network: Input BiologicalNetwork (reference or empirical).
        target_genes: Optional subset of genes to evaluate (defaults to all genes in network).
        gene_symbol_map: Optional mapping of Ensembl ID -> gene symbol.

    Returns:
        Dictionary mapping canonical Ensembl gene ID -> GeneNetworkTopologyMetrics.
    """
    symbols = gene_symbol_map or {}
    G = build_networkx_graph(network)

    if target_genes is not None:
        eval_genes = [g.strip().upper() for g in target_genes if g]
    else:
        eval_genes = list(G.nodes())

    all_degrees = dict(G.degree())
    all_deg_values = list(all_degrees.values()) if all_degrees else [0]

    # Precompute clustering coefficients
    clustering_dict = nx.clustering(G) if G.number_of_nodes() > 0 else {}

    results: dict[str, GeneNetworkTopologyMetrics] = {}

    for gene in eval_genes:
        sym = symbols.get(gene)
        if gene in G:
            deg = int(all_degrees.get(gene, 0))
            clust = float(clustering_dict.get(gene, 0.0))
            neighbors = list(G.neighbors(gene))
            n_size = len(neighbors)
            is_iso = (deg == 0)
            pct = float(stats_percentile(all_deg_values, deg))
        else:
            # Gene is not in the network (e.g. unannotated in reference network)
            deg = 0
            clust = 0.0
            n_size = 0
            is_iso = True
            pct = 0.0

        metrics = GeneNetworkTopologyMetrics(
            canonical_ensembl_gene_id=gene,
            gene_symbol=sym,
            network_id=network.network_id,
            degree=deg,
            clustering_coefficient=clust,
            neighborhood_size=n_size,
            is_isolated=is_iso,
            degree_percentile=pct,
            reference_degree_bias_note=REFERENCE_DEGREE_BIAS_DISCLAIMER,
        )
        results[gene] = metrics

    return results


def stats_percentile(all_values: Sequence[int], value: int) -> float:
    """Calculate empirical percentile rank of a value within a list of values."""
    if not all_values:
        return 0.0
    arr = np.array(all_values)
    return float(np.mean(arr <= value) * 100.0)


def extract_chemical_responsive_subnetwork(
    network: BiologicalNetwork,
    responsive_genes: Iterable[str],
    subnetwork_id: str | None = None,
    subnetwork_name: str | None = None,
) -> BiologicalNetwork:
    """Extract induced subnetwork connecting chemical-responsive genes.

    Preserves exact edge metadata from the parent network while isolating
    interactions among observed differentially responsive genes.
    """
    resp_set = {g.strip().upper() for g in responsive_genes if g}
    sub_edges = tuple(
        e for e in network.edges
        if e.gene_a_canonical_id.upper() in resp_set and e.gene_b_canonical_id.upper() in resp_set
    )

    s_id = subnetwork_id or f"{network.network_id}_induced_subnetwork"
    s_name = subnetwork_name or f"{network.network_name} (Induced Subnetwork: {len(resp_set)} genes)"

    prov = dict(network.provenance)
    prov["induced_subnetwork_parent"] = network.network_id
    prov["responsive_genes_count"] = len(resp_set)
    prov["subnetwork_edges_count"] = len(sub_edges)

    return BiologicalNetwork(
        network_id=s_id,
        network_name=s_name,
        relationship_type=network.relationship_type,
        edges=sub_edges,
        species=network.species,
        taxon_id=network.taxon_id,
        source_database=network.source_database,
        database_version=network.database_version,
        qc_status=network.qc_status,
        scientific_limitations=(
            f"Induced subnetwork connecting {len(resp_set)} chemical-responsive genes. "
            f"Edge connectivity does NOT imply direct chemical target interaction."
        ),
        provenance=prov,
    )


def detect_network_modules(
    network: BiologicalNetwork,
    min_module_size: int = 3,
) -> dict[str, list[str]]:
    """Detect discrete functional modules/communities in a network using greedy modularity.

    Args:
        network: Input BiologicalNetwork.
        min_module_size: Minimum number of genes required to retain a module.

    Returns:
        Dictionary mapping module ID (e.g. 'module_1') to list of canonical gene IDs.
    """
    G = build_networkx_graph(network)
    if G.number_of_edges() == 0:
        return {}

    # Run greedy modularity community detection
    try:
        communities = nx.community.greedy_modularity_communities(G)
    except Exception:
        # Fallback to connected components
        communities = list(nx.connected_components(G))

    modules: dict[str, list[str]] = {}
    mod_idx = 1
    for comm in communities:
        gene_list = sorted(list(comm))
        if len(gene_list) >= min_module_size:
            modules[f"module_{mod_idx}"] = gene_list
            mod_idx += 1

    return modules
