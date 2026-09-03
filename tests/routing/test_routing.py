"""Format detection and routing acceptance tests."""

import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from src.ingestion import AssetStatus, Route
from src.routing import route_metadata_assets

from .fixtures import add_asset, routing_metadata, write_asset


class RoutingTests(unittest.TestCase):
    def route(self, path: Path, declared: str, modality: str = "bulk_rna_seq"):
        return route_metadata_assets(routing_metadata(path, declared, modality))[0]

    def test_valid_fastq(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.fastq", "fastq"), "FASTQ")
        self.assertEqual((AssetStatus.READY, Route.RAW_READ_PIPELINE, "FASTQ"), (result.status, result.route, result.detected_format))

    def test_valid_fq(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.fq", "fastq"), "FASTQ")
        self.assertEqual(AssetStatus.READY, result.status)

    def test_valid_gzip_fastq(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.fastq.gz", "fastq_gz"), "FASTQ")
        self.assertEqual(AssetStatus.READY, result.status)
        self.assertIn("gzip compression", result.evidence)

    def test_malformed_fastq_is_invalid(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "bad.fastq", "malformed_fastq"), "FASTQ")
        self.assertEqual(AssetStatus.INVALID, result.status)
        self.assertIn("malformed_fastq", {item.rule for item in result.findings})

    def test_empty_fastq_is_invalid(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "empty.fastq", "empty"), "FASTQ")
        self.assertEqual(AssetStatus.INVALID, result.status)

    def test_sra_signature_routes_raw_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.sra", "sra"), "SRA")
        self.assertEqual((AssetStatus.READY, Route.RAW_READ_PIPELINE), (result.status, result.route))

    def test_sra_extension_without_signature_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.sra", "unknown"), "SRA")
        self.assertEqual(AssetStatus.NEEDS_REVIEW, result.status)

    def test_integer_raw_counts_with_bulk_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "integer_matrix"), "raw_counts")
        self.assertEqual((AssetStatus.READY, Route.BULK_RAW_COUNT_PIPELINE), (result.status, result.route))

    def test_non_integer_matrix_is_processed_expression(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "processed_matrix"), "processed_expression_matrix")
        self.assertEqual((AssetStatus.READY, Route.PROCESSED_EXPRESSION_PIPELINE), (result.status, result.route))

    def test_generic_integer_matrix_is_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "integer_matrix"), "other")
        self.assertEqual((AssetStatus.NEEDS_REVIEW, Route.MANUAL_REVIEW), (result.status, result.route))

    def test_umi_counts_with_single_cell_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "integer_matrix"), "UMI_counts", "scRNA_seq")
        self.assertEqual((AssetStatus.READY, Route.SINGLE_CELL_COUNT_PIPELINE), (result.status, result.route))

    def test_umi_counts_with_bulk_metadata_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "integer_matrix"), "UMI_counts")
        self.assertEqual((AssetStatus.NEEDS_REVIEW, Route.MANUAL_REVIEW), (result.status, result.route))

    def test_h5ad_recognition(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "cells.h5ad", "h5ad"), "H5AD", "scRNA_seq")
        self.assertEqual("H5AD", result.detected_format)

    def test_h5ad_with_single_cell_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "cells.h5ad", "h5ad"), "H5AD", "scRNA_seq")
        self.assertEqual((AssetStatus.READY, Route.SINGLE_CELL_OBJECT_PIPELINE), (result.status, result.route))

    def test_h5ad_with_bulk_metadata_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "cells.h5ad", "h5ad"), "H5AD")
        self.assertEqual((AssetStatus.NEEDS_REVIEW, Route.MANUAL_REVIEW), (result.status, result.route))

    def test_rds_is_conservative(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "object.rds", "rds"), "RDS", "scRNA_seq")
        self.assertEqual(("RDS", AssetStatus.NEEDS_REVIEW, Route.MANUAL_REVIEW), (result.detected_format, result.status, result.route))

    def test_unknown_file_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "asset.bin", "unknown"), "other")
        self.assertEqual(AssetStatus.NEEDS_REVIEW, result.status)

    def test_missing_file_is_invalid(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(Path(directory) / "missing.fastq", "FASTQ")
        self.assertEqual(AssetStatus.INVALID, result.status)

    def test_empty_matrix_is_invalid_when_zero_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "empty"), "raw_counts")
        self.assertEqual(AssetStatus.INVALID, result.status)

    def test_header_only_matrix_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "empty_matrix"), "raw_counts")
        self.assertEqual(AssetStatus.NEEDS_REVIEW, result.status)

    def test_declared_and_detected_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "reads.fastq", "fastq"), "FASTQ")
        self.assertEqual(result.declared_format, result.detected_format)
        self.assertNotIn("declared_detected_mismatch", {item.rule for item in result.findings})

    def test_declared_detected_conflict_needs_review(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "matrix.tsv", "processed_matrix"), "raw_counts")
        self.assertEqual((AssetStatus.NEEDS_REVIEW, Route.MANUAL_REVIEW), (result.status, result.route))
        self.assertIn("declared_detected_mismatch", {item.rule for item in result.findings})

    def test_multiple_assets_for_one_assay(self):
        with tempfile.TemporaryDirectory() as directory:
            first = write_asset(Path(directory) / "reads_1.fastq", "fastq")
            second = write_asset(Path(directory) / "reads_2.fastq", "fastq")
            document = add_asset(routing_metadata(first, "FASTQ"), second, "asset_2")
            results = route_metadata_assets(document)
        self.assertEqual(2, len(results))
        self.assertTrue(all(item.status == AssetStatus.READY for item in results))

    def test_repeated_results_are_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            document = routing_metadata(write_asset(Path(directory) / "matrix.tsv", "integer_matrix"), "other")
            first = route_metadata_assets(document)
            second = route_metadata_assets(document)
        self.assertEqual(first, second)

    def test_input_file_is_not_modified(self):
        with tempfile.TemporaryDirectory() as directory:
            path = write_asset(Path(directory) / "reads.fastq", "fastq")
            before = path.read_bytes()
            self.route(path, "FASTQ")
            after = path.read_bytes()
        self.assertEqual(before, after)

    def test_canonical_metadata_is_not_modified(self):
        with tempfile.TemporaryDirectory() as directory:
            document = routing_metadata(write_asset(Path(directory) / "reads.fastq", "fastq"), "FASTQ")
            before = deepcopy(document)
            route_metadata_assets(document)
        self.assertEqual(before, document)

    def test_filename_words_do_not_determine_format(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.route(write_asset(Path(directory) / "raw_counts_expression.bin", "unknown"), "other")
        self.assertIsNone(result.detected_format)


if __name__ == "__main__":
    unittest.main()
