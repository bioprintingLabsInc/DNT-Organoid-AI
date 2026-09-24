"""Unit tests for gene identifier type resolution and non-guessing safeguards."""

import unittest

from src.gene_harmonization.config import load_config
from src.gene_harmonization.models import GeneIdentifierType, Severity
from src.gene_harmonization.type_detector import (
    normalize_declared_type,
    resolve_identifier_type,
)


class TypeDetectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = load_config()
        cls.supported = cls.config["supported_identifier_types"]

    def test_alias_normalization(self):
        self.assertEqual(
            normalize_declared_type("Ensembl Gene ID", self.supported),
            GeneIdentifierType.ENSEMBL_GENE_ID,
        )
        self.assertEqual(
            normalize_declared_type("ENSG", self.supported),
            GeneIdentifierType.ENSEMBL_GENE_ID,
        )
        self.assertEqual(
            normalize_declared_type("Gene Symbol", self.supported),
            GeneIdentifierType.GENE_SYMBOL,
        )
        self.assertEqual(
            normalize_declared_type("symbol", self.supported),
            GeneIdentifierType.GENE_SYMBOL,
        )
        self.assertEqual(
            normalize_declared_type("Entrez Gene ID", self.supported),
            GeneIdentifierType.ENTREZ_GENE_ID,
        )
        self.assertEqual(
            normalize_declared_type("HGNC ID", self.supported),
            GeneIdentifierType.HGNC_ID,
        )

    def test_unsupported_declared_type_yields_error(self):
        gene_ids = ("G1", "G2")
        id_type, findings = resolve_identifier_type(
            gene_ids, "unsupported_foreign_id", self.supported, "asset_1"
        )
        self.assertIsNone(id_type)
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, Severity.ERROR)
        self.assertEqual(findings[0].rule_id, "unsupported_identifier_type")

    def test_declared_type_mismatch_warning(self):
        # Declared Ensembl but given gene symbols
        gene_ids = ("TP53", "GAPDH", "ACTB")
        id_type, findings = resolve_identifier_type(
            gene_ids, "Ensembl Gene ID", self.supported, "asset_1"
        )
        self.assertEqual(id_type, GeneIdentifierType.ENSEMBL_GENE_ID)
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, Severity.REVIEW)
        self.assertEqual(findings[0].rule_id, "identifier_type_mismatch")

    def test_undeclared_unambiguous_ensembl_detection(self):
        gene_ids = ("ENSG00000141510", "ENSG00000111640.8", "ENSG00000075624")
        id_type, findings = resolve_identifier_type(
            gene_ids, None, self.supported, "asset_1"
        )
        self.assertEqual(id_type, GeneIdentifierType.ENSEMBL_GENE_ID)
        self.assertEqual(len(findings), 0)

    def test_undeclared_pure_numbers_refuses_to_guess(self):
        # Could be Entrez or row indices: must not silently guess
        gene_ids = ("7157", "2597", "60")
        id_type, findings = resolve_identifier_type(
            gene_ids, None, self.supported, "asset_1"
        )
        self.assertIsNone(id_type)
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, Severity.REVIEW)
        self.assertEqual(findings[0].rule_id, "ambiguous_identifier_type")

    def test_empty_gene_ids_yields_error(self):
        gene_ids = ("", "   ")
        id_type, findings = resolve_identifier_type(
            gene_ids, "Ensembl Gene ID", self.supported, "asset_1"
        )
        self.assertIsNone(id_type)
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].severity, Severity.ERROR)
        self.assertEqual(findings[0].rule_id, "empty_gene_identifiers")


if __name__ == "__main__":
    unittest.main()
