"""Unit tests for gene reference catalog loading, indexing, and identifier lookups."""

import unittest
from pathlib import Path

from src.gene_harmonization.config import load_config
from src.gene_harmonization.models import GeneIdentifierType, GeneMappingStatus
from src.gene_harmonization.reference import GeneReference


class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = load_config()
        cls.reference = GeneReference.load_from_config(cls.config)

    def test_reference_metadata_and_checksum(self):
        self.assertEqual(self.reference.reference_id, "ensembl_human_genes")
        self.assertEqual(self.reference.reference_version, "112")
        self.assertEqual(self.reference.genome_build, "GRCh38.p14")
        self.assertEqual(len(self.reference.checksum), 64)
        self.assertGreater(len(self.reference.records), 50)

    def test_plain_ensembl_lookup(self):
        res = self.reference.lookup_ensembl("ENSG00000141510")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(res.approved_symbol, "TP53")
        self.assertIsNone(res.reason)

    def test_versioned_ensembl_lookup(self):
        res = self.reference.lookup_ensembl("ENSG00000141510.16")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(res.approved_symbol, "TP53")

    def test_unmapped_ensembl_lookup(self):
        res = self.reference.lookup_ensembl("ENSG99999999999")
        self.assertEqual(res.status, GeneMappingStatus.UNMAPPED)
        self.assertIsNone(res.canonical_gene_id)
        self.assertIsNotNone(res.reason)
        self.assertIn("not found in reference", res.reason)

    def test_ensembl_gene_without_hgnc_annotation(self):
        # ENSG00000228037 is a canonical Ensembl 112 gene with no HGNC cross-reference
        res = self.reference.lookup_ensembl("ENSG00000228037")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000228037")
        self.assertIsNone(res.approved_symbol)

    def test_approved_symbol_lookup(self):
        res = self.reference.lookup_symbol("GAPDH")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000111640")
        self.assertEqual(res.approved_symbol, "GAPDH")

    def test_symbol_case_insensitivity(self):
        res = self.reference.lookup_symbol("gapdh")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000111640")
        self.assertEqual(res.approved_symbol, "GAPDH")

    def test_synonym_lookup(self):
        # TBR2 is an established synonym for EOMES
        res = self.reference.lookup_symbol("TBR2")
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000163508")
        self.assertEqual(res.approved_symbol, "EOMES")

    def test_ambiguous_symbol_lookup_preserves_candidates(self):
        # FLIP is a multi-hit synonym matching multiple canonical Ensembl genes
        res = self.reference.lookup_symbol("FLIP")
        self.assertEqual(res.status, GeneMappingStatus.AMBIGUOUS)
        self.assertIsNone(res.canonical_gene_id)
        self.assertIsNone(res.approved_symbol)
        self.assertEqual(len(res.candidate_canonical_ids), 2)
        self.assertIn("ENSG00000000460", res.candidate_canonical_ids)
        self.assertIn("ENSG00000003402", res.candidate_canonical_ids)
        self.assertIn("matches multiple canonical Ensembl IDs", res.reason)

    def test_ambiguous_synonym_lookup(self):
        # SRC2 is a multi-hit synonym matching multiple canonical Ensembl genes
        res = self.reference.lookup_symbol("SRC2")
        self.assertEqual(res.status, GeneMappingStatus.AMBIGUOUS)
        self.assertIsNone(res.canonical_gene_id)
        self.assertIn("ENSG00000000938", res.candidate_canonical_ids)
        self.assertIn("ENSG00000140396", res.candidate_canonical_ids)

    def test_unmapped_symbol_lookup(self):
        res = self.reference.lookup_symbol("UNKNOWN_SYMBOL_12345")
        self.assertEqual(res.status, GeneMappingStatus.UNMAPPED)
        self.assertIsNone(res.canonical_gene_id)
        self.assertIn("not found in reference", res.reason)

    def test_empty_or_whitespace_lookup_is_invalid(self):
        res = self.reference.lookup("   ", GeneIdentifierType.GENE_SYMBOL)
        self.assertEqual(res.status, GeneMappingStatus.INVALID)
        self.assertIsNone(res.canonical_gene_id)

    def test_entrez_id_lookup(self):
        # 7157 -> TP53 (ENSG00000141510)
        res = self.reference.lookup("7157", GeneIdentifierType.ENTREZ_GENE_ID)
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(res.approved_symbol, "TP53")

    def test_hgnc_id_lookup(self):
        # HGNC:11998 -> TP53 (ENSG00000141510)
        res = self.reference.lookup("HGNC:11998", GeneIdentifierType.HGNC_ID)
        self.assertEqual(res.status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(res.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(res.approved_symbol, "TP53")


if __name__ == "__main__":
    unittest.main()
