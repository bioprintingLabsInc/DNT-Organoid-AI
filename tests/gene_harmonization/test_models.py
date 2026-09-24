"""Unit tests for gene harmonization data models, immutability, and status semantics."""

import unittest
from src.gene_harmonization.models import (
    CanonicalCollision,
    Finding,
    GeneIdentifierType,
    GeneMappingStatus,
    HarmonizationSummary,
    HarmonizedDataset,
    HarmonizedGene,
    Severity,
    Status,
    ordered,
    status_for,
)
from src.qc.matrix import CountMatrix


class ModelsTests(unittest.TestCase):
    def test_gene_mapping_status_values(self):
        self.assertEqual(GeneMappingStatus.UNIQUELY_MAPPED, "UNIQUELY_MAPPED")
        self.assertEqual(GeneMappingStatus.AMBIGUOUS, "AMBIGUOUS")
        self.assertEqual(GeneMappingStatus.UNMAPPED, "UNMAPPED")
        self.assertEqual(GeneMappingStatus.INVALID, "INVALID")

    def test_harmonized_gene_immutability(self):
        gene = HarmonizedGene(
            source_index=0,
            original_gene_id="ENSG00000141510.15",
            identifier_type="ensembl_gene_id",
            canonical_gene_id="ENSG00000141510",
            approved_symbol="TP53",
            mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
            mapping_reason=None,
            candidate_canonical_ids=("ENSG00000141510",),
            reference_id="ensembl_human_genes",
            reference_version="112",
            source_asset_id="asset_1",
        )
        with self.assertRaises((TypeError, AttributeError)):
            gene.canonical_gene_id = "ENSG99999999999"  # type: ignore

    def test_canonical_collision_immutability(self):
        col = CanonicalCollision(
            canonical_gene_id="ENSG00000141510",
            original_gene_ids=("ENSG00000141510.15", "ENSG00000141510.16"),
            source_indices=(0, 1),
            approved_symbol="TP53",
        )
        self.assertEqual(len(col.original_gene_ids), 2)
        with self.assertRaises((TypeError, AttributeError)):
            col.canonical_gene_id = "ENSG00000111640"  # type: ignore

    def test_finding_ordering_and_precedence(self):
        f_info = Finding(Severity.INFO, "info_rule", "gene", "G1", "p", "info msg")
        f_warn = Finding(Severity.WARNING, "warn_rule", "gene", "G2", "p", "warn msg")
        f_rev = Finding(Severity.REVIEW, "rev_rule", "gene", "G3", "p", "rev msg")
        f_err = Finding(Severity.ERROR, "err_rule", "gene", "G4", "p", "err msg")

        raw_list = [f_info, f_err, f_warn, f_rev]
        sorted_findings = ordered(raw_list)

        self.assertEqual(sorted_findings[0].severity, Severity.ERROR)
        self.assertEqual(sorted_findings[1].severity, Severity.REVIEW)
        self.assertEqual(sorted_findings[2].severity, Severity.WARNING)
        self.assertEqual(sorted_findings[3].severity, Severity.INFO)

    def test_status_for_severity_precedence(self):
        f_info = Finding(Severity.INFO, "r", "e", "1", "p", "m")
        f_warn = Finding(Severity.WARNING, "r", "e", "2", "p", "m")
        f_rev = Finding(Severity.REVIEW, "r", "e", "3", "p", "m")
        f_err = Finding(Severity.ERROR, "r", "e", "4", "p", "m")

        self.assertEqual(status_for([]), Status.PASS)
        self.assertEqual(status_for([f_info]), Status.PASS)
        self.assertEqual(status_for([f_info, f_warn]), Status.PASS_WITH_WARNINGS)
        self.assertEqual(status_for([f_warn, f_rev]), Status.NEEDS_REVIEW)
        self.assertEqual(status_for([f_rev, f_err]), Status.FAIL)

    def test_harmonized_dataset_to_dict_and_slots(self):
        summary = HarmonizationSummary(1, 1, 0, 0, 0, 0, 0, 100.0, 0.0)
        matrix = CountMatrix(gene_ids=("G1",), sample_ids=("S1",), columns=((10,),))
        dataset = HarmonizedDataset(
            status=Status.PASS,
            findings=(),
            summary=summary,
            genes=(),
            collisions=(),
            matrix=matrix,
            source_asset_id="asset_1",
            assay_id="assay_1",
            sample_ids=("S1",),
            reference_identity="ensembl_human_genes",
            reference_version="112",
            reference_checksum="sha256:abc",
            config_version="0.1.0",
            config_checksum="sha256:def",
            harmonization_version="bulk_gene_harmonization_v1",
        )
        d = dataset.to_dict()
        self.assertEqual(d["status"], "PASS")
        self.assertEqual(d["source_asset_id"], "asset_1")
        self.assertEqual(d["reference_identity"], "ensembl_human_genes")
        self.assertEqual(d["reference_version"], "112")


if __name__ == "__main__":
    unittest.main()
