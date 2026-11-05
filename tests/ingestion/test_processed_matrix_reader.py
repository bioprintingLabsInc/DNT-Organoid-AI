"""Tests for Phase U3.2 Processed Expression Matrix Intake."""

import gzip
from pathlib import Path
import tempfile
import unittest

from src.ingestion.processed_matrix_reader import (
    ProcessedExpressionMatrix,
    ProcessedMatrixData,
    read_processed_matrix,
)


class ProcessedMatrixReaderTests(unittest.TestCase):
    """Test reading, structural validation, and numerical preservation of processed expression matrices."""

    def test_read_valid_tsv_matrix(self):
        """Read a valid TSV processed expression matrix preserving float values and headers."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene_id\tsample1\tsample2\n")
            f.write("ENSG00000141510\t12.54321\t0.0\n")
            f.write("ENSG00000012048\t105.78901\t24.11005\n")
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertTrue(data.is_valid)
            self.assertIsNotNone(data.matrix)
            m = data.matrix

            self.assertEqual(2, m.number_of_genes)
            self.assertEqual(2, m.number_of_samples)
            self.assertEqual(("ENSG00000141510", "ENSG00000012048"), m.gene_ids)
            self.assertEqual(("sample1", "sample2"), m.sample_ids)

            # Exact floating point preservation (no rounding or integer truncation)
            self.assertEqual(12.54321, m.get_value("ENSG00000141510", "sample1"))
            self.assertEqual(0.0, m.get_value("ENSG00000141510", "sample2"))
            self.assertEqual(105.78901, m.get_value("ENSG00000012048", "sample1"))
            self.assertEqual(24.11005, m.get_value("ENSG00000012048", "sample2"))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_read_valid_csv_gz_matrix(self):
        """Read a gzip-compressed CSV processed expression matrix."""
        with tempfile.NamedTemporaryFile("wb", suffix=".csv.gz", delete=False) as f:
            content = "gene_id,sample_a,sample_b\nTP53,5.5,10.25\nEGFR,0.12,0.48\n"
            f.write(gzip.compress(content.encode("utf-8")))
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertTrue(data.is_valid)
            self.assertIsNotNone(data.matrix)
            self.assertEqual(2, data.matrix.number_of_genes)
            self.assertEqual(2, data.matrix.number_of_samples)
            self.assertEqual(5.5, data.matrix.get_value("TP53", "sample_a"))
            self.assertEqual(0.48, data.matrix.get_value("EGFR", "sample_b"))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_empty_file_rejected(self):
        """Empty file emits ERROR finding and returns invalid data."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("")
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertFalse(data.is_valid)
            self.assertIsNone(data.matrix)
            self.assertTrue(any(f.rule == "empty_file" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_nonexistent_file_rejected(self):
        """Nonexistent file emits ERROR finding."""
        data = read_processed_matrix(Path("nonexistent_matrix_file_12345.tsv"))
        self.assertFalse(data.is_valid)
        self.assertIsNone(data.matrix)
        self.assertTrue(any(f.rule == "file_not_found" for f in data.findings))

    def test_duplicate_sample_identifiers_rejected(self):
        """Duplicate sample columns in header emit unique_sample_identifiers_required ERROR."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\tGSM1\tGSM1\n")
            f.write("GENE1\t1.0\t2.0\n")
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertFalse(data.is_valid)
            self.assertIsNone(data.matrix)
            self.assertTrue(any(f.rule == "unique_sample_identifiers_required" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_malformed_row_width_rejected(self):
        """Row with missing or extra columns emits malformed_row ERROR."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\tsample1\tsample2\n")
            f.write("GENE1\t1.0\t2.0\n")
            f.write("GENE2\t3.0\n")  # Only 2 columns instead of 3
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertFalse(data.is_valid)
            self.assertIsNone(data.matrix)
            self.assertTrue(any(f.rule == "malformed_row" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_missing_values_rejected(self):
        """Missing or NA values emit missing_expression_values_forbidden ERROR."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\tsample1\tsample2\n")
            f.write("GENE1\t1.0\t\n")  # Empty cell
            f.write("GENE2\tNA\t2.5\n")  # NA string
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertFalse(data.is_valid)
            self.assertIsNone(data.matrix)
            self.assertTrue(any(f.rule == "missing_expression_values_forbidden" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_non_numeric_values_rejected(self):
        """Non-numeric string values emit numeric_expression_required ERROR."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\tsample1\tsample2\n")
            f.write("GENE1\t1.0\tnot_a_number\n")
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertFalse(data.is_valid)
            self.assertIsNone(data.matrix)
            self.assertTrue(any(f.rule == "numeric_expression_required" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_duplicate_gene_identifiers_preserved_with_review_finding(self):
        """Duplicate gene rows are preserved in source row order without silent deduplication."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\ts1\ts2\n")
            f.write("DUPLICATE_GENE\t1.5\t2.5\n")
            f.write("UNIQUE_GENE\t3.0\t4.0\n")
            f.write("DUPLICATE_GENE\t5.5\t6.5\n")
            tmp_path = Path(f.name)

        try:
            data = read_processed_matrix(tmp_path)
            self.assertTrue(data.is_valid)
            self.assertIsNotNone(data.matrix)
            m = data.matrix

            # Both duplicate rows are preserved
            self.assertEqual(3, m.number_of_genes)
            self.assertEqual(("DUPLICATE_GENE", "UNIQUE_GENE", "DUPLICATE_GENE"), m.gene_ids)
            self.assertEqual((1.5, 3.0, 5.5), m.columns[0])
            self.assertEqual((2.5, 4.0, 6.5), m.columns[1])

            # Emits REVIEW finding for duplicate gene
            self.assertTrue(any(f.rule == "duplicate_gene_identifiers" and f.severity == "REVIEW" for f in data.findings))
        finally:
            tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
