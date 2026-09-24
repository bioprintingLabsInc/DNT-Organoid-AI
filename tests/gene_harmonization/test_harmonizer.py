"""Comprehensive integration and unit tests for Gene Harmonization v1."""

import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from src.gene_harmonization import (
    GeneMappingStatus,
    Status,
    harmonize_bulk_counts,
)
from src.qc import Status as QCStatus, validate_bulk_counts
from src.qc.matrix import CountMatrix

try:
    from .fixtures import harmonization_metadata, write_test_matrix
except ImportError:
    from tests.gene_harmonization.fixtures import harmonization_metadata, write_test_matrix


class HarmonizerTests(unittest.TestCase):
    def test_clean_pass_with_plain_ensembl_ids(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.PASS)
        self.assertEqual(result.summary.total_genes, 3)
        self.assertEqual(result.summary.uniquely_mapped_genes, 3)
        self.assertEqual(result.summary.unmapped_genes, 0)
        self.assertEqual(result.summary.ambiguous_genes, 0)
        self.assertEqual(result.summary.canonical_collision_genes, 0)
        self.assertEqual(len(result.genes), 3)

        # First gene TP53
        g1 = result.genes[0]
        self.assertEqual(g1.original_gene_id, "ENSG00000141510")
        self.assertEqual(g1.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(g1.approved_symbol, "TP53")
        self.assertEqual(g1.mapping_status, GeneMappingStatus.UNIQUELY_MAPPED)
        self.assertEqual(g1.source_index, 0)

    def test_clean_pass_with_versioned_ensembl_ids(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "ENSG00000141510.16\t10\t20\t30\t40\t50\t60\n"
                "ENSG00000111640.8\t100\t110\t120\t130\t140\t150\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.PASS)
        self.assertEqual(result.summary.uniquely_mapped_genes, 2)
        # Verify version suffix stripped for canonical ID, but original identifier string preserved
        self.assertEqual(result.genes[0].original_gene_id, "ENSG00000141510.16")
        self.assertEqual(result.genes[0].canonical_gene_id, "ENSG00000141510")
        self.assertEqual(result.genes[0].approved_symbol, "TP53")
        self.assertEqual(result.genes[1].original_gene_id, "ENSG00000111640.8")
        self.assertEqual(result.genes[1].canonical_gene_id, "ENSG00000111640")
        self.assertEqual(result.genes[1].approved_symbol, "GAPDH")

    def test_clean_pass_with_valid_gene_symbols(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "TP53\t10\t20\t30\t40\t50\t60\n"
                "GAPDH\t100\t110\t120\t130\t140\t150\n"
                "SOX2\t50\t60\t70\t80\t90\t100\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="gene_symbol")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.PASS)
        self.assertEqual(result.summary.uniquely_mapped_genes, 3)
        self.assertEqual(result.genes[0].canonical_gene_id, "ENSG00000141510")
        self.assertEqual(result.genes[1].canonical_gene_id, "ENSG00000111640")
        self.assertEqual(result.genes[2].canonical_gene_id, "ENSG00000181449")
        self.assertEqual(result.genes[2].approved_symbol, "SOX2")

    def test_unmapped_identifiers_retained_and_surface_review(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "TP53\t10\t20\t30\t40\t50\t60\n"
                "NONEXISTENT_GENE_XYZ\t5\t5\t5\t5\t5\t5\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="gene_symbol")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        # Without threshold, unmapped gene is retained with INFO finding
        self.assertEqual(result.summary.total_genes, 2)
        self.assertEqual(result.summary.uniquely_mapped_genes, 1)
        self.assertEqual(result.summary.unmapped_genes, 1)
        self.assertEqual(result.summary.percentage_uniquely_mapped, 50.0)
        self.assertEqual(result.summary.percentage_unresolved, 50.0)
        self.assertEqual(result.genes[1].original_gene_id, "NONEXISTENT_GENE_XYZ")
        self.assertEqual(result.genes[1].mapping_status, GeneMappingStatus.UNMAPPED)
        self.assertIsNone(result.genes[1].canonical_gene_id)
        self.assertIn("not found in reference", result.genes[1].mapping_reason or "")
        self.assertIn("unmapped_gene_identifier", {f.rule_id for f in result.findings})

    def test_ambiguous_identifiers_preserve_candidates_and_never_force_unique(self):
        with tempfile.TemporaryDirectory() as d:
            # FLIP is a multi-hit synonym with multiple candidate Ensembl IDs
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "TP53\t10\t20\t30\t40\t50\t60\n"
                "FLIP\t2\t2\t2\t2\t2\t2\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="gene_symbol")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.NEEDS_REVIEW)
        self.assertEqual(result.summary.ambiguous_genes, 1)
        flip_gene = result.genes[1]
        self.assertEqual(flip_gene.mapping_status, GeneMappingStatus.AMBIGUOUS)
        self.assertIsNone(flip_gene.canonical_gene_id)
        self.assertEqual(len(flip_gene.candidate_canonical_ids), 2)
        self.assertIn("ENSG00000000460", flip_gene.candidate_canonical_ids)
        self.assertIn("ENSG00000003402", flip_gene.candidate_canonical_ids)
        self.assertIn("ambiguous_gene_identifier", {f.rule_id for f in result.findings})

    def test_multiple_source_ids_mapping_to_one_canonical_gene_collision(self):
        with tempfile.TemporaryDirectory() as d:
            # Two versioned Ensembl IDs that resolve to the same canonical Ensembl ID
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "ENSG00000141510.15\t10\t20\t30\t40\t50\t60\n"
                "ENSG00000141510.16\t15\t25\t35\t45\t55\t65\n"
                "ENSG00000111640.8\t100\t110\t120\t130\t140\t150\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.NEEDS_REVIEW)
        self.assertEqual(result.summary.total_genes, 3)
        self.assertEqual(result.summary.canonical_collision_genes, 1)
        self.assertEqual(result.summary.colliding_source_gene_count, 2)
        self.assertEqual(len(result.collisions), 1)

        collision = result.collisions[0]
        self.assertEqual(collision.canonical_gene_id, "ENSG00000141510")
        self.assertEqual(collision.original_gene_ids, ("ENSG00000141510.15", "ENSG00000141510.16"))
        self.assertEqual(collision.source_indices, (0, 1))

        # Check that NO aggregation was performed: 3 separate rows exist
        self.assertEqual(len(result.matrix.columns[0]), 3)
        self.assertEqual(result.matrix.columns[0][0], 10)
        self.assertEqual(result.matrix.columns[0][1], 15)
        self.assertIn("canonical_gene_collision", {f.rule_id for f in result.findings})

    def test_duplicate_source_identifiers_surface_review_and_collision(self):
        with tempfile.TemporaryDirectory() as d:
            content = (
                "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
                "TP53\t10\t20\t30\t40\t50\t60\n"
                "TP53\t99\t88\t77\t66\t55\t44\n"
            )
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="gene_symbol")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(result.status, Status.NEEDS_REVIEW)
        self.assertEqual(result.summary.total_genes, 2)
        self.assertEqual(result.summary.canonical_collision_genes, 1)
        self.assertEqual(result.summary.colliding_source_gene_count, 2)

        rule_ids = {f.rule_id for f in result.findings}
        self.assertIn("duplicate_source_gene_identifiers", rule_ids)
        self.assertIn("canonical_gene_collision", rule_ids)

        # Verify neither row was dropped or summed
        self.assertEqual(len(result.matrix.columns[0]), 2)
        self.assertEqual(result.matrix.columns[0][0], 10)
        self.assertEqual(result.matrix.columns[0][1], 99)

    def test_no_expression_values_are_lost_or_altered(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            result = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        expected_cols = (
            (10, 100, 500),
            (20, 110, 520),
            (30, 120, 510),
            (40, 130, 530),
            (50, 140, 540),
            (60, 150, 550),
        )
        self.assertEqual(result.matrix.columns, expected_cols)

    def test_deterministic_output_and_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            res_a = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")
            res_b = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(res_a, res_b)
        self.assertEqual(res_a.reference_identity, "ensembl_human_genes")
        self.assertEqual(res_a.reference_version, "112")
        self.assertEqual(len(res_a.reference_checksum), 64)
        self.assertEqual(len(res_a.config_checksum), 64)
        self.assertEqual(res_a.harmonization_version, "bulk_gene_harmonization_v1")

    def test_source_matrix_and_metadata_remain_immutable(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            before_file = matrix_path.read_bytes()
            before_meta = deepcopy(meta)

            harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

            self.assertEqual(matrix_path.read_bytes(), before_file)
            self.assertEqual(meta, before_meta)

    def test_qc_gate_failure_stops_harmonization(self):
        with tempfile.TemporaryDirectory() as d:
            # Create a count matrix with negative values that fails QC
            content = "gene\tsample_c1\nG1\t-1\n"
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            qc_res = validate_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")
            self.assertEqual(qc_res.status, QCStatus.FAIL)

            harm_res = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", qc_result=qc_res
            )

        self.assertEqual(harm_res.status, Status.FAIL)
        self.assertIn("qc_status_failure", {f.rule_id for f in harm_res.findings})
        self.assertIsNone(harm_res.matrix)

    def test_qc_gate_needs_review_unauthorized_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.qc.result import DatasetResult
            qc_res = DatasetResult(
                status=QCStatus.NEEDS_REVIEW,
                findings=(),
                number_of_genes=3,
                number_of_samples=1,
                source_asset_id="counts_asset",
                assay_id="bulk_assay",
                study_ids=(),
                experiment_ids=(),
                ruleset_version="1.0",
                ruleset_checksum="abc",
                calculation_version="1.0",
                samples=(),
                conditions=(),
            )
            harm_res = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", qc_result=qc_res
            )

        self.assertEqual(harm_res.status, Status.FAIL)
        self.assertIn("unauthorized_upstream_qc_needs_review", {f.rule_id for f in harm_res.findings})
        self.assertIsNone(harm_res.matrix)

    def test_qc_gate_bare_boolean_fails(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.qc.result import DatasetResult
            qc_res = DatasetResult(
                status=QCStatus.NEEDS_REVIEW,
                findings=(),
                number_of_genes=3,
                number_of_samples=1,
                source_asset_id="counts_asset",
                assay_id="bulk_assay",
                study_ids=(),
                experiment_ids=(),
                ruleset_version="1.0",
                ruleset_checksum="abc",
                calculation_version="1.0",
                samples=(),
                conditions=(),
            )
            # Passing bare allow_needs_review_qc=True without disposition MUST fail
            harm_res = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", qc_result=qc_res, allow_needs_review_qc=True
            )

        self.assertEqual(harm_res.status, Status.FAIL)
        self.assertIn("unauthorized_upstream_qc_needs_review", {f.rule_id for f in harm_res.findings})
        self.assertIsNone(harm_res.matrix)

    def test_qc_gate_structured_disposition_object_proceeds(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.qc.result import DatasetResult
            qc_res = DatasetResult(
                status=QCStatus.NEEDS_REVIEW,
                findings=(),
                number_of_genes=3,
                number_of_samples=1,
                source_asset_id="counts_asset",
                assay_id="bulk_assay",
                study_ids=(),
                experiment_ids=(),
                ruleset_version="1.0",
                ruleset_checksum="abc",
                calculation_version="1.0",
                samples=(),
                conditions=(),
            )
            from src.gene_harmonization import QCReviewDisposition
            disp = QCReviewDisposition(
                decision="ACCEPTED",
                reason="Minor sample correlation warning manually inspected and cleared by PI",
                reviewer="Scientist-42",
                review_date="2026-09-24T12:00:00Z",
                affected_asset_id="counts_asset",
            )
            harm_res = harmonize_bulk_counts(
                matrix_path,
                meta,
                "counts_asset",
                "bulk_assay",
                qc_result=qc_res,
                qc_disposition=disp,
            )

        self.assertEqual(harm_res.status, Status.NEEDS_REVIEW)
        self.assertIn("upstream_qc_needs_review_authorized", {f.rule_id for f in harm_res.findings})
        self.assertIsNotNone(harm_res.matrix)
        self.assertIsNotNone(harm_res.qc_review_disposition)
        self.assertEqual(harm_res.qc_review_disposition.decision, "ACCEPTED")
        self.assertEqual(harm_res.qc_review_disposition.reviewer, "Scientist-42")

    def test_ensembl_gene_without_hgnc_annotation_in_matrix(self):
        with tempfile.TemporaryDirectory() as d:
            content = "gene\ts1\ts2\nENSG00000228037\t10\t20\n"
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            harm_res = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(harm_res.status, Status.PASS)
        self.assertEqual(harm_res.summary.uniquely_mapped_genes, 1)
        self.assertEqual(harm_res.summary.unmapped_genes, 0)
        gene = harm_res.genes[0]
        self.assertEqual(gene.canonical_gene_id, "ENSG00000228037")
        self.assertIsNone(gene.approved_symbol)
        self.assertEqual(gene.mapping_status, GeneMappingStatus.UNIQUELY_MAPPED)

    def test_qc_gate_needs_review_authorized_by_disposition_proceeds(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.qc.result import DatasetResult
            qc_res = DatasetResult(
                status=QCStatus.NEEDS_REVIEW,
                findings=(),
                number_of_genes=3,
                number_of_samples=1,
                source_asset_id="counts_asset",
                assay_id="bulk_assay",
                study_ids=(),
                experiment_ids=(),
                ruleset_version="1.0",
                ruleset_checksum="abc",
                calculation_version="1.0",
                samples=(),
                conditions=(),
            )
            harm_res = harmonize_bulk_counts(
                matrix_path,
                meta,
                "counts_asset",
                "bulk_assay",
                qc_result=qc_res,
                qc_disposition="REVIEWED_AND_CONFIRMED_VIABLE",
            )

        self.assertEqual(harm_res.status, Status.NEEDS_REVIEW)
        self.assertIn("upstream_qc_needs_review_authorized", {f.rule_id for f in harm_res.findings})
        self.assertIsNotNone(harm_res.matrix)
        self.assertEqual(harm_res.upstream_qc_status, "NEEDS_REVIEW")

    def test_qc_gate_pass_with_warnings_proceeds(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.qc.result import DatasetResult
            qc_res = DatasetResult(
                status=QCStatus.PASS_WITH_WARNINGS,
                findings=(),
                number_of_genes=3,
                number_of_samples=1,
                source_asset_id="counts_asset",
                assay_id="bulk_assay",
                study_ids=(),
                experiment_ids=(),
                ruleset_version="1.0",
                ruleset_checksum="abc",
                calculation_version="1.0",
                samples=(),
                conditions=(),
            )
            harm_res = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", qc_result=qc_res
            )

        self.assertEqual(harm_res.status, Status.PASS_WITH_WARNINGS)
        self.assertIn("upstream_qc_pass_with_warnings", {f.rule_id for f in harm_res.findings})
        self.assertIsNotNone(harm_res.matrix)
        self.assertEqual(harm_res.upstream_qc_status, "PASS_WITH_WARNINGS")

    def test_coverage_threshold_and_percentages(self):
        with tempfile.TemporaryDirectory() as d:
            # 2 genes: 1 mapped, 1 unmapped -> 50% uniquely mapped
            content = "gene\tsample1\nENSG00000141510\t10\nUNKNOWN_GENE\t20\n"
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv", content)
            meta = harmonization_metadata(matrix_path, gene_identifier_type="Ensembl Gene ID")

            from src.gene_harmonization.config import load_config
            cfg = load_config()
            cfg_review = deepcopy(cfg)
            cfg_review["coverage_thresholds"]["minimum_mapping_coverage_percentage"] = 80.0

            harm_res = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", config=cfg_review
            )

            self.assertEqual(harm_res.status, Status.NEEDS_REVIEW)
            self.assertIn("low_mapping_coverage", {f.rule_id for f in harm_res.findings})
            self.assertEqual(harm_res.summary.total_genes, 2)
            self.assertEqual(harm_res.summary.uniquely_mapped_genes, 1)
            self.assertEqual(harm_res.summary.percentage_uniquely_mapped, 50.0)
            self.assertEqual(harm_res.summary.percentage_unresolved, 50.0)

            # Test configured ERROR severity
            cfg_error = deepcopy(cfg_review)
            cfg_error["structural_rules"]["low_mapping_coverage"] = "ERROR"
            harm_res_err = harmonize_bulk_counts(
                matrix_path, meta, "counts_asset", "bulk_assay", config=cfg_error
            )
            self.assertEqual(harm_res_err.status, Status.FAIL)

    def test_unsupported_identifier_type_fails(self):
        with tempfile.TemporaryDirectory() as d:
            matrix_path = write_test_matrix(Path(d) / "matrix.tsv")
            meta = harmonization_metadata(matrix_path, gene_identifier_type="unknown_unsupported_format")

            harm_res = harmonize_bulk_counts(matrix_path, meta, "counts_asset", "bulk_assay")

        self.assertEqual(harm_res.status, Status.FAIL)
        self.assertIn("unsupported_identifier_type", {f.rule_id for f in harm_res.findings})

    def test_accepts_preloaded_count_matrix(self):
        matrix = CountMatrix(
            gene_ids=("ENSG00000141510", "ENSG00000111640"),
            sample_ids=("S1", "S2"),
            columns=((10, 20), (30, 40)),
        )
        meta = harmonization_metadata(Path("virtual.tsv"), gene_identifier_type="Ensembl Gene ID")

        result = harmonize_bulk_counts(matrix, meta, "counts_asset", "bulk_assay")
        self.assertEqual(result.status, Status.PASS)
        self.assertEqual(result.summary.uniquely_mapped_genes, 2)
        self.assertEqual(result.matrix, matrix)


if __name__ == "__main__":
    unittest.main()
