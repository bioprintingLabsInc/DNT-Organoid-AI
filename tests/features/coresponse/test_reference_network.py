"""Unit and integration tests for biological reference network ingestion in Stage 8B."""

from pathlib import Path
import pandas as pd
import pytest

from src.features.coresponse.models import (
    CoresponseQCStatus,
    EvidenceCategory,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
)
from src.features.coresponse.reference_network import (
    filter_network_by_confidence,
    filter_network_by_genes,
    import_reference_network_from_table,
)


def test_import_reference_network_from_table(tmp_path: Path) -> None:
    """Test importing reference network preserving interaction nature, confidence, and version."""
    tsv_content = (
        "gene1\tgene2\tsym1\tsym2\tnature\tevidence\tscore\n"
        "ENSG00000141510\tENSG00000146648\tTP53\tEGFR\tPHYSICAL_BINDING\tEXPERIMENTAL\t0.95\n"
        "ENSG00000141510\tENSG00000105383\tTP53\tCDKN1A\tREGULATORY\tDATABASE_CURATED\t0.88\n"
        "ENSG00000146648\tENSG00000157764\tEGFR\tBRAF\tFUNCTIONAL\tCOMBINED\t0.75\n"
        "ENSG00000157764\tENSG00000171862\tBRAF\tPTEN\tPREDICTED\tTEXTMINING\t0.42\n"
    )
    tsv_file = tmp_path / "reference_ppi_v12.tsv"
    tsv_file.write_text(tsv_content, encoding="utf-8")

    net = import_reference_network_from_table(
        table_path=tsv_file,
        source_database="STRING",
        database_version="v12.0",
        gene_a_col="gene1",
        gene_b_col="gene2",
        gene_a_symbol_col="sym1",
        gene_b_symbol_col="sym2",
        interaction_nature_col="nature",
        evidence_category_col="evidence",
        confidence_score_col="score",
    )

    assert net.species == "Homo sapiens"
    assert net.taxon_id == 9606
    assert net.source_database == "STRING"
    assert net.database_version == "v12.0"
    assert net.edge_count == 4
    assert len(net.unique_genes) == 5

    # Check TP53 - EGFR physical interaction
    tp53_edges = net.get_edges_for_gene("ENSG00000141510")
    assert len(tp53_edges) == 2

    egfr_edge = [e for e in tp53_edges if "ENSG00000146648" in (e.gene_a_canonical_id, e.gene_b_canonical_id)][0]
    assert egfr_edge.interaction_nature == InteractionNature.PHYSICAL_INTERACTION
    assert egfr_edge.is_direct_binding is True
    assert egfr_edge.confidence_score == 0.95
    assert egfr_edge.relationship_type == NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION
    assert egfr_edge.partition_assignment == PartitionAssignment.UNRESTRICTED_REFERENCE

    # CRITICAL INVARIANT: Reference interaction is NEVER marked as chemical-specific
    assert egfr_edge.is_chemical_specific is False
    assert egfr_edge.causality_claim == "REFERENCE_BIOLOGICAL_ASSOCIATION"


def test_reference_network_multi_version_support_without_code_changes(tmp_path: Path) -> None:
    """Test importing different versions of reference databases without code modifications."""
    df_v1 = pd.DataFrame({
        "gene_a": ["ENSG001", "ENSG002"],
        "gene_b": ["ENSG002", "ENSG003"],
        "score": [0.8, 0.9],
    })
    net_v1 = import_reference_network_from_table(
        table_path=df_v1,
        source_database="BioGRID",
        database_version="4.4.230",
        confidence_score_col="score",
    )
    assert net_v1.database_version == "4.4.230"
    assert net_v1.edge_count == 2

    # Import newer version
    df_v2 = pd.DataFrame({
        "gene_a": ["ENSG001", "ENSG002", "ENSG003"],
        "gene_b": ["ENSG002", "ENSG004", "ENSG004"],
        "score": [0.85, 0.92, 0.77],
    })
    net_v2 = import_reference_network_from_table(
        table_path=df_v2,
        source_database="BioGRID",
        database_version="4.4.235",
        confidence_score_col="score",
    )
    assert net_v2.database_version == "4.4.235"
    assert net_v2.edge_count == 3


def test_reference_network_confidence_filtering() -> None:
    """Test filtering network by minimum confidence threshold."""
    df = pd.DataFrame({
        "gene_a": ["ENSG001", "ENSG001", "ENSG002"],
        "gene_b": ["ENSG002", "ENSG003", "ENSG003"],
        "score": [950, 450, 700],  # 0-1000 scale
    })
    net = import_reference_network_from_table(
        table_path=df,
        source_database="STRING",
        database_version="v12.0",
        confidence_score_col="score",
    )
    assert net.edge_count == 3
    # Scores scaled to [0, 1]
    assert net.edges[0].confidence_score == 0.95

    # Filter with min_confidence = 0.70
    filtered = filter_network_by_confidence(net, min_confidence=0.70)
    assert filtered.edge_count == 2
    assert all(e.confidence_score >= 0.70 for e in filtered.edges)


def test_genes_absent_from_reference_network() -> None:
    """Test that genes absent from reference network are handled cleanly without fabricating edges."""
    df = pd.DataFrame({
        "gene_a": ["ENSG001"],
        "gene_b": ["ENSG002"],
    })
    net = import_reference_network_from_table(
        table_path=df,
        source_database="TestDB",
        database_version="1.0",
    )
    # Query edges for gene not in reference
    absent_edges = net.get_edges_for_gene("ENSG999999")
    assert len(absent_edges) == 0

    # Subnetwork filtering with absent genes
    sub = filter_network_by_genes(net, target_genes=["ENSG999999", "ENSG888888"])
    assert sub.edge_count == 0
