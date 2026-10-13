"""Unit tests for Stage 8C EnsemblToEntrezAdapter mapping."""

from pathlib import Path
import pytest

from src.differential_expression.models import GeneDifferentialExpressionResult
from src.features.pathway.mapping import EnsemblToEntrezAdapter
from src.features.pathway.models import MappingStatus
from src.gene_harmonization.reference import GeneReference, ReferenceRecord


@pytest.fixture
def mock_reference() -> GeneReference:
    """Fixture providing a controlled in-memory GeneReference."""
    records = (
        ReferenceRecord(
            ensembl_gene_id="ENSG00000000003",
            gene_symbol="TSPAN6",
            entrez_gene_id="7105",
            hgnc_id="HGNC:11858",
            gene_biotype="protein_coding",
            synonyms=("T245", "TM4SF6"),
            chromosome="X",
            description="tetraspanin 6",
        ),
        ReferenceRecord(
            ensembl_gene_id="ENSG00000000005",
            gene_symbol="TNMD",
            entrez_gene_id="",  # Unmapped
            hgnc_id="HGNC:17757",
            gene_biotype="protein_coding",
            synonyms=(),
            chromosome="X",
            description="tenomodulin",
        ),
        ReferenceRecord(
            ensembl_gene_id="ENSG00000240758",
            gene_symbol="HILPDA-AS1",
            entrez_gene_id="125775236",  # Duplicate Entrez 1
            hgnc_id="HGNC:54321",
            gene_biotype="lncRNA",
            synonyms=(),
            chromosome="1",
            description="HILPDA antisense RNA 1",
        ),
        ReferenceRecord(
            ensembl_gene_id="ENSG00000272601",
            gene_symbol="EFCAB3P1",
            entrez_gene_id="125775236",  # Duplicate Entrez 2
            hgnc_id="HGNC:54322",
            gene_biotype="pseudogene",
            synonyms=(),
            chromosome="1",
            description="EF-hand calcium binding domain 3 pseudogene 1",
        ),
        ReferenceRecord(
            ensembl_gene_id="ENSG00000999999",
            gene_symbol="MULTI",
            entrez_gene_id="999|888",  # One-to-many
            hgnc_id="HGNC:99999",
            gene_biotype="protein_coding",
            synonyms=(),
            chromosome="2",
            description="multi gene",
        ),
    )
    return GeneReference(
        reference_id="test_ref_v1",
        reference_version="Test_v1.0",
        genome_build="GRCh38",
        source_authority="Ensembl",
        records=records,
        checksum="test_checksum_1234",
    )


def test_valid_ensembl_to_entrez_mapping(mock_reference: GeneReference) -> None:
    adapter = EnsemblToEntrezAdapter(gene_reference=mock_reference)
    genes = [
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000003",
            approved_symbol="TSPAN6",
            base_mean=100.0,
            log2_fold_change=1.5,
            lfc_standard_error=0.3,
            wald_statistic=5.0,
            p_value=1e-6,
            adjusted_p_value_bh=1e-5,
            result_status="OK",
        )
    ]
    res = adapter.map_genes(genes)
    assert res.total_input_genes == 1
    assert res.mapped_count == 1
    assert res.unmapped_count == 0
    rec = res.mappings[0]
    assert rec.canonical_gene_id == "ENSG00000000003"
    assert rec.entrez_gene_id == "7105"
    assert rec.approved_symbol == "TSPAN6"
    assert rec.mapping_status == MappingStatus.MAPPED.value


def test_unmapped_gene_handling(mock_reference: GeneReference) -> None:
    adapter = EnsemblToEntrezAdapter(gene_reference=mock_reference)
    genes = [
        # Gene with empty Entrez in reference
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000005",
            approved_symbol="TNMD",
            base_mean=50.0,
            log2_fold_change=-0.5,
            lfc_standard_error=0.2,
            wald_statistic=-2.5,
            p_value=0.01,
            adjusted_p_value_bh=0.05,
            result_status="OK",
        ),
        # Gene completely missing from reference
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00009999999",
            approved_symbol="UNKNOWN",
            base_mean=0.0,
            log2_fold_change=0.0,
            lfc_standard_error=1.0,
            wald_statistic=0.0,
            p_value=1.0,
            adjusted_p_value_bh=1.0,
            result_status="OK",
        ),
    ]
    res = adapter.map_genes(genes)
    assert res.unmapped_count == 2
    assert res.mappings[0].mapping_status == MappingStatus.UNMAPPED.value
    assert res.mappings[0].entrez_gene_id is None
    assert res.mappings[1].mapping_status == MappingStatus.UNMAPPED.value
    assert res.mappings[1].entrez_gene_id is None


def test_duplicate_and_one_to_many_mapping(mock_reference: GeneReference) -> None:
    adapter = EnsemblToEntrezAdapter(gene_reference=mock_reference)
    genes = [
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000240758",
            approved_symbol="HILPDA-AS1",
            base_mean=20.0,
            log2_fold_change=2.0,
            lfc_standard_error=0.5,
            wald_statistic=4.0,
            p_value=1e-4,
            adjusted_p_value_bh=1e-3,
            result_status="OK",
        ),
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000272601",
            approved_symbol="EFCAB3P1",
            base_mean=15.0,
            log2_fold_change=1.0,
            lfc_standard_error=0.5,
            wald_statistic=2.0,
            p_value=0.04,
            adjusted_p_value_bh=0.1,
            result_status="OK",
        ),
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000999999",
            approved_symbol="MULTI",
            base_mean=30.0,
            log2_fold_change=-1.0,
            lfc_standard_error=0.2,
            wald_statistic=-5.0,
            p_value=1e-6,
            adjusted_p_value_bh=1e-5,
            result_status="OK",
        ),
    ]
    res = adapter.map_genes(genes)
    assert res.duplicate_entrez_count == 2
    assert res.one_to_many_count == 1
    assert res.mappings[0].mapping_status == MappingStatus.DUPLICATE_ENTREZ.value
    assert res.mappings[1].mapping_status == MappingStatus.DUPLICATE_ENTREZ.value
    assert res.mappings[2].mapping_status == MappingStatus.ONE_TO_MANY.value
    assert res.mappings[2].entrez_gene_id == "999"


def test_deterministic_mapping_and_provenance(mock_reference: GeneReference) -> None:
    adapter = EnsemblToEntrezAdapter(gene_reference=mock_reference)
    genes = [
        GeneDifferentialExpressionResult(
            canonical_gene_id="ENSG00000000003",
            approved_symbol="TSPAN6",
            base_mean=10.0,
            log2_fold_change=1.0,
            lfc_standard_error=0.5,
            wald_statistic=2.0,
            p_value=0.05,
            adjusted_p_value_bh=0.1,
            result_status="OK",
        )
    ]
    res1 = adapter.map_genes(genes)
    res2 = adapter.map_genes(genes)
    assert res1.to_dict() == res2.to_dict()
    assert res1.reference_version == "Test_v1.0"
    assert res1.reference_checksum == "test_checksum_1234"


def test_real_ensembl_reference_loads_and_maps() -> None:
    adapter = EnsemblToEntrezAdapter()
    assert "112" in adapter.gene_reference.reference_version
    assert "ENSG00000105851" in adapter.gene_reference._by_ensembl
    rec = adapter.gene_reference._by_ensembl["ENSG00000105851"]
    assert rec.gene_symbol == "PIK3CG"
    assert rec.entrez_gene_id == "5294"
