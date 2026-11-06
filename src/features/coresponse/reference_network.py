"""Stage 8B: Import and Ingestion of Established Biological Reference Networks.

Supports universal, versioned import of human gene and protein interaction networks
(e.g., STRING, BioGRID, IntAct, or generic tabular interactome resources) while:
1. Enforcing canonical Ensembl gene identities.
2. Preserving interaction type, physical vs. functional nature, evidence category, and confidence.
3. Preserving species (Homo sapiens, taxon 9606) and explicit reference version.
4. Enforcing the invariant that reference interactions are NOT chemical-specific co-responses.
"""

from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import pandas as pd

from src.features.coresponse.models import (
    BiologicalNetwork,
    CoresponseQCStatus,
    EvidenceCategory,
    GeneNetworkEdge,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
    derive_edge_id,
)


def import_reference_network_from_table(
    table_path: str | Path | pd.DataFrame,
    source_database: str,
    database_version: str,
    gene_a_col: str = "gene_a",
    gene_b_col: str = "gene_b",
    gene_a_symbol_col: str | None = None,
    gene_b_symbol_col: str | None = None,
    interaction_nature_col: str | None = None,
    evidence_category_col: str | None = None,
    confidence_score_col: str | None = None,
    default_interaction_nature: InteractionNature = InteractionNature.FUNCTIONAL_ASSOCIATION,
    default_evidence_category: EvidenceCategory = EvidenceCategory.DATABASE_CURATED,
    network_id: str | None = None,
    network_name: str | None = None,
    min_confidence: float | None = None,
    species: str = "Homo sapiens",
    taxon_id: int = 9606,
    provenance_metadata: dict[str, Any] | None = None,
) -> BiologicalNetwork:
    """Universal tabular importer for established biological reference networks.

    Args:
        table_path: File path (CSV/TSV) or DataFrame containing interaction records.
        source_database: Name of external database (e.g. 'STRING', 'BioGRID', 'IntAct').
        database_version: Version identifier of the external resource (e.g. 'v12.0', '4.4.234').
        gene_a_col: Column containing Ensembl ID for interactor A.
        gene_b_col: Column containing Ensembl ID for interactor B.
        gene_a_symbol_col: Optional column for approved gene symbol for interactor A.
        gene_b_symbol_col: Optional column for approved gene symbol for interactor B.
        interaction_nature_col: Optional column containing interaction nature string.
        evidence_category_col: Optional column containing evidence category string.
        confidence_score_col: Optional column containing numeric confidence score.
        default_interaction_nature: Default nature when column is absent or value is unmapped.
        default_evidence_category: Default category when column is absent or value is unmapped.
        network_id: Optional network identifier.
        network_name: Optional human-readable network title.
        min_confidence: Optional confidence cutoff threshold.
        species: Organism species string (defaults to 'Homo sapiens').
        taxon_id: NCBI Taxonomy ID (defaults to 9606).
        provenance_metadata: Optional additional metadata dict.

    Returns:
        Immutable BiologicalNetwork instance with typed GeneNetworkEdge records.
    """
    if isinstance(table_path, pd.DataFrame):
        df = table_path.copy()
    else:
        path = Path(table_path)
        if not path.exists():
            raise FileNotFoundError(f"Reference network file not found at '{path}'.")
        sep = "\t" if path.suffix in {".tsv", ".tab", ".txt"} else ","
        df = pd.read_csv(path, sep=sep, low_memory=False)

    if gene_a_col not in df.columns or gene_b_col not in df.columns:
        raise ValueError(
            f"Required gene columns '{gene_a_col}' or '{gene_b_col}' not found in reference table. "
            f"Available columns: {list(df.columns)}"
        )

    edges: list[GeneNetworkEdge] = []
    seen_pairs: set[tuple[str, str]] = set()

    for _, row in df.iterrows():
        g1_raw = str(row[gene_a_col]).strip() if pd.notna(row[gene_a_col]) else ""
        g2_raw = str(row[gene_b_col]).strip() if pd.notna(row[gene_b_col]) else ""

        if not g1_raw or not g2_raw or g1_raw.upper() == g2_raw.upper():
            # Skip empty or self-interacting loop edges
            continue

        # Canonical ordering for undirected network
        g1, g2 = sorted([g1_raw, g2_raw])
        pair_key = (g1.upper(), g2.upper())
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        sym1 = str(row[gene_a_symbol_col]).strip() if gene_a_symbol_col and pd.notna(row.get(gene_a_symbol_col)) else None
        sym2 = str(row[gene_b_symbol_col]).strip() if gene_b_symbol_col and pd.notna(row.get(gene_b_symbol_col)) else None

        # Resolve confidence score
        conf = None
        if confidence_score_col and pd.notna(row.get(confidence_score_col)):
            try:
                conf = float(row[confidence_score_col])
                # Scale 0-1000 integer scores (e.g. STRING) to [0.0, 1.0] if > 1.0
                if conf > 1.0 and conf <= 1000.0:
                    conf = conf / 1000.0
            except (ValueError, TypeError):
                conf = None

        if min_confidence is not None and conf is not None and conf < min_confidence:
            continue

        # Resolve Interaction Nature
        nature = default_interaction_nature
        if interaction_nature_col and pd.notna(row.get(interaction_nature_col)):
            val = str(row[interaction_nature_col]).strip().upper()
            if "PHYSICAL" in val or "BINDING" in val or "COMPLEX" in val:
                nature = InteractionNature.PHYSICAL_INTERACTION
            elif "PREDICT" in val or "HOMOLOGY" in val or "TEXTMINING" in val:
                nature = InteractionNature.PREDICTED_INTERACTION
            elif "REGULATORY" in val or "TRANSCRIPTION" in val:
                nature = InteractionNature.REGULATORY_INTERACTION
            elif "FUNCTIONAL" in val:
                nature = InteractionNature.FUNCTIONAL_ASSOCIATION

        # Resolve Evidence Category
        cat = default_evidence_category
        if evidence_category_col and pd.notna(row.get(evidence_category_col)):
            val = str(row[evidence_category_col]).strip().upper()
            if "EXP" in val or "ASSAY" in val or "Y2H" in val or "MS" in val:
                cat = EvidenceCategory.EXPERIMENTAL_ASSAY
            elif "CURAT" in val or "DATABASE" in val:
                cat = EvidenceCategory.DATABASE_CURATED
            elif "TEXT" in val:
                cat = EvidenceCategory.TEXT_MINING
            elif "COMP" in val or "PRED" in val:
                cat = EvidenceCategory.COMPUTATIONAL_PREDICTION
            elif "COMBINED" in val:
                cat = EvidenceCategory.COMBINED_SCORE

        is_phys = (nature == InteractionNature.PHYSICAL_INTERACTION)
        edge_id = derive_edge_id(g1, g2, NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION, source_database)

        edge = GeneNetworkEdge(
            edge_id=edge_id,
            gene_a_canonical_id=g1,
            gene_b_canonical_id=g2,
            gene_a_symbol=sym1,
            gene_b_symbol=sym2,
            relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
            interaction_nature=nature,
            evidence_source=source_database,
            evidence_version=database_version,
            evidence_category=cat,
            confidence_score=conf,
            partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
            qc_status=CoresponseQCStatus.VALID,
            is_direct_binding=is_phys,
            is_chemical_specific=False,  # CRITICAL: Reference interaction is NOT chemical-specific co-response
            causality_claim="REFERENCE_BIOLOGICAL_ASSOCIATION",
            scientific_limitations=(
                "Established reference network association. Reflects known prior biological or physical interaction. "
                "Does NOT establish that these genes co-respond to chemical exposure or share directional regulatory causality."
            ),
            provenance_notes=f"Source: {source_database} {database_version}; Taxon: {taxon_id} ({species})",
        )
        edges.append(edge)

    net_id = network_id or f"ref_{source_database.lower()}_{database_version.lower()}"
    net_name = network_name or f"{source_database} Reference Interactome ({database_version})"

    prov = {
        "source_database": source_database,
        "database_version": database_version,
        "species": species,
        "taxon_id": taxon_id,
        "total_edges": len(edges),
        "min_confidence_filter": min_confidence,
        "imported_at": datetime.now(timezone.utc).isoformat(),
    }
    if provenance_metadata:
        prov.update(provenance_metadata)

    return BiologicalNetwork(
        network_id=net_id,
        network_name=net_name,
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        edges=tuple(edges),
        species=species,
        taxon_id=taxon_id,
        source_database=source_database,
        database_version=database_version,
        qc_status=CoresponseQCStatus.VALID,
        scientific_limitations="Independent biological reference network. Unbiased by empirical exposure cohorts.",
        provenance=prov,
    )


def filter_network_by_confidence(network: BiologicalNetwork, min_confidence: float) -> BiologicalNetwork:
    """Filter reference network edges by minimum confidence threshold."""
    filtered_edges = tuple(
        e for e in network.edges
        if e.confidence_score is not None and e.confidence_score >= min_confidence
    )
    new_prov = dict(network.provenance)
    new_prov["confidence_threshold_applied"] = min_confidence
    new_prov["filtered_edges_count"] = len(filtered_edges)

    return BiologicalNetwork(
        network_id=f"{network.network_id}_conf_{min_confidence}",
        network_name=f"{network.network_name} (score >= {min_confidence})",
        relationship_type=network.relationship_type,
        edges=filtered_edges,
        species=network.species,
        taxon_id=network.taxon_id,
        source_database=network.source_database,
        database_version=network.database_version,
        qc_status=network.qc_status,
        scientific_limitations=network.scientific_limitations,
        provenance=new_prov,
    )


def filter_network_by_genes(network: BiologicalNetwork, target_genes: Iterable[str]) -> BiologicalNetwork:
    """Extract induced subnetwork containing only edges where both genes are in target_genes."""
    target_set = {g.strip().upper() for g in target_genes if g}
    sub_edges = tuple(
        e for e in network.edges
        if e.gene_a_canonical_id.upper() in target_set and e.gene_b_canonical_id.upper() in target_set
    )
    new_prov = dict(network.provenance)
    new_prov["gene_filter_target_count"] = len(target_set)
    new_prov["subnetwork_edges_count"] = len(sub_edges)

    return BiologicalNetwork(
        network_id=f"{network.network_id}_subnetwork",
        network_name=f"{network.network_name} (Subnetwork: {len(target_set)} genes)",
        relationship_type=network.relationship_type,
        edges=sub_edges,
        species=network.species,
        taxon_id=network.taxon_id,
        source_database=network.source_database,
        database_version=network.database_version,
        qc_status=network.qc_status,
        scientific_limitations=network.scientific_limitations,
        provenance=new_prov,
    )
