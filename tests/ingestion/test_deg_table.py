"""Unit and edge case tests for published differential-expression summary table ingestion (Phase U3.1)."""

import gzip
from pathlib import Path
import tempfile
import unittest

from src.ingestion.deg_table_reader import (
    DegTableData,
    discover_deg_columns,
    read_deg_table,
)
from tests.routing.test_phase_u1 import create_test_xlsx


class DegTableReaderTests(unittest.TestCase):
    """Verify safe, bounded, and non-fabricating published DEG table parsing."""

    def test_read_valid_tsv_deg_table(self):
        """Clean TSV with standard DESeq2 headers is parsed with exact numeric values."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "results.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\tgene_name\n"
                "ENSG00000101017\t3.397\t2.62e-106\t7.21e-102\tCD40\n"
                "ENSG00000130876\t-1.850\t0.00045\t0.0032\tSLC7A10\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertEqual(2, data.total_rows)
            self.assertEqual("ENSG00000101017", data.rows[0].raw_gene_id)
            self.assertEqual("CD40", data.rows[0].raw_gene_symbol)
            self.assertAlmostEqual(3.397, data.rows[0].log2_fold_change)
            self.assertAlmostEqual(2.62e-106, data.rows[0].p_value)
            self.assertAlmostEqual(7.21e-102, data.rows[0].adjusted_p_value)

            # Confirm negative fold change preserved exactly
            self.assertAlmostEqual(-1.850, data.rows[1].log2_fold_change)

    def test_read_valid_csv_deg_table(self):
        """Comma-delimited DEG table with limma-style headers is detected and parsed."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "limma_results.csv"
            p.write_text(
                "symbol,logFC,AveExpr,t,P.Value,adj.P.Val\n"
                "GAPDH,2.5,12.4,5.8,1.2e-5,4.5e-4\n"
                "ACTB,-0.8,14.1,-2.3,0.021,0.089\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertEqual(2, data.total_rows)
            self.assertEqual("GAPDH", data.rows[0].raw_gene_id)
            self.assertAlmostEqual(2.5, data.rows[0].log2_fold_change)
            self.assertAlmostEqual(12.4, data.rows[0].base_mean)
            self.assertAlmostEqual(5.8, data.rows[0].wald_statistic)
            self.assertAlmostEqual(1.2e-5, data.rows[0].p_value)
            self.assertAlmostEqual(4.5e-4, data.rows[0].adjusted_p_value)

    def test_read_gzip_deg_table(self):
        """Gzip-compressed TSV DEG table is uncompressed in-memory and parsed."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "results.tsv.gz"
            content = "gene_id\tlog2FoldChange\tpvalue\tpadj\nENSG1\t1.5\t0.01\t0.05\n"
            with gzip.open(p, "wt", encoding="utf-8") as f:
                f.write(content)

            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertEqual(1, data.total_rows)
            self.assertEqual("ENSG1", data.rows[0].raw_gene_id)
            self.assertTrue(data.provenance["is_gzip"])

    def test_read_xlsx_deg_table_single_sheet(self):
        """Excel workbook with a single DEG sheet is detected and parsed without external libraries."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "published_supp.xlsx"
            create_test_xlsx(p, {
                "VPAvsCTRL": [
                    ["gene_id", "log2FoldChange", "pvalue", "padj", "gene_name"],
                    ["ENSG00000101017", 3.397, 2.62e-106, 7.21e-102, "CD40"],
                    ["ENSG00000130876", -1.85, 0.00045, 0.0032, "SLC7A10"],
                ],
                "Readme": [
                    ["Title", "Study summary"],
                    ["Organoid", "Dorsal forebrain"],
                ],
            })

            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertEqual(2, data.total_rows)
            self.assertEqual("VPAvsCTRL", data.sheet_name)
            self.assertEqual("ENSG00000101017", data.rows[0].raw_gene_id)

    def test_read_xlsx_deg_table_explicit_sheet_name(self):
        """Caller can specify explicit sheet_name to target specific contrast sheet."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "multi_contrast.xlsx"
            create_test_xlsx(p, {
                "Dose_10uM": [
                    ["gene_id", "log2FoldChange", "padj"],
                    ["G1", 1.2, 0.01],
                ],
                "Dose_50uM": [
                    ["gene_id", "log2FoldChange", "padj"],
                    ["G1", 2.8, 0.0001],
                ],
            })

            data = read_deg_table(p, sheet_name="Dose_50uM")
            self.assertTrue(data.is_valid)
            self.assertEqual("Dose_50uM", data.sheet_name)
            self.assertAlmostEqual(2.8, data.rows[0].log2_fold_change)

    def test_ambiguous_xlsx_sheets_requires_explicit_selection(self):
        """Workbook with multiple DEG candidate sheets fails without explicit sheet_name."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "ambiguous_deg.xlsx"
            create_test_xlsx(p, {
                "Dose_10uM": [["gene_id", "log2FoldChange", "padj"], ["G1", 1.2, 0.01]],
                "Dose_50uM": [["gene_id", "log2FoldChange", "padj"], ["G1", 2.8, 0.0001]],
            })

            data = read_deg_table(p)
            self.assertFalse(data.is_valid)
            self.assertTrue(any(f.rule == "ambiguous_deg_sheets" for f in data.findings))

    def test_zero_data_fabrication_missing_se_wald_basemean_remain_none(self):
        """CRITICAL: Missing statistics are NEVER invented; they remain strictly None."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "minimal.tsv"
            p.write_text(
                "gene\tlog2fc\tpadj\n"
                "A1BG\t1.4\t0.02\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            row = data.rows[0]
            self.assertIsNone(row.lfc_standard_error)
            self.assertIsNone(row.wald_statistic)
            self.assertIsNone(row.base_mean)
            self.assertIsNone(row.p_value)  # Only padj was published

    def test_na_and_nan_pvalues_preserved_as_none_not_dropped(self):
        """Genes with NA/NaN p-values (DESeq2 outlier filtering) are preserved as None, not dropped."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "outliers.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "ENSG1\t3.5\tNA\tNA\n"
                "ENSG2\t-2.1\t0.001\t0.02\n"
                "ENSG3\t0.0\tNaN\tNaN\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertEqual(3, data.total_rows)

            self.assertAlmostEqual(3.5, data.rows[0].log2_fold_change)
            self.assertIsNone(data.rows[0].p_value)
            self.assertIsNone(data.rows[0].adjusted_p_value)

            self.assertAlmostEqual(-2.1, data.rows[1].log2_fold_change)
            self.assertAlmostEqual(0.001, data.rows[1].p_value)

            self.assertAlmostEqual(0.0, data.rows[2].log2_fold_change)
            self.assertIsNone(data.rows[2].p_value)

    def test_invalid_pvalue_bounds_flagged_and_set_none(self):
        """P-values outside [0, 1] range trigger warning finding and are set to None."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "invalid_p.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "G1\t1.5\t1.5\t0.05\n"
                "G2\t-1.0\t0.01\t-0.2\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertTrue(data.is_valid)
            self.assertTrue(any(f.rule == "invalid_significance_values" for f in data.findings))
            self.assertIsNone(data.rows[0].p_value)
            self.assertIsNone(data.rows[1].adjusted_p_value)

    def test_missing_gene_column_rejected(self):
        """Table missing any gene identifier column is rejected with missing_gene_column error."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "no_gene.tsv"
            p.write_text("log2FoldChange\tpvalue\tpadj\n1.5\t0.01\t0.05\n", encoding="utf-8")
            data = read_deg_table(p)
            self.assertFalse(data.is_valid)
            self.assertTrue(any(f.rule == "missing_gene_column" for f in data.findings))

    def test_missing_effect_size_rejected(self):
        """Table missing effect size column is rejected with missing_effect_size_column error."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "no_effect.tsv"
            p.write_text("gene_id\tpvalue\tpadj\nENSG1\t0.01\t0.05\n", encoding="utf-8")
            data = read_deg_table(p)
            self.assertFalse(data.is_valid)
            self.assertTrue(any(f.rule == "missing_effect_size_column" for f in data.findings))

    def test_explicit_column_mapping_overrides(self):
        """Caller can provide custom column mapping for non-standard header formats."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "custom_headers.tsv"
            p.write_text(
                "CustomIdentifier\tDifference\tSignificance\n"
                "MY_GENE\t2.4\t0.0005\n",
                encoding="utf-8",
            )
            data = read_deg_table(
                p,
                column_mapping={"gene_id": "CustomIdentifier", "log2fc": "Difference", "padj": "Significance"},
            )
            self.assertTrue(data.is_valid)
            self.assertEqual("MY_GENE", data.rows[0].raw_gene_id)
            self.assertAlmostEqual(2.4, data.rows[0].log2_fold_change)
            self.assertAlmostEqual(0.0005, data.rows[0].adjusted_p_value)


    def test_ordinary_fold_change_rejected_without_explicit_confirmation(self):
        """Ambiguous linear fold change (e.g. fold_change, fc) is NEVER automatically treated as log2FC."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "linear_fc.tsv"
            p.write_text(
                "gene_id\tfold_change\tpvalue\tpadj\n"
                "ENSG1\t2.0\t0.01\t0.05\n",
                encoding="utf-8",
            )
            data = read_deg_table(p)
            self.assertFalse(data.is_valid)
            self.assertTrue(any(f.rule == "ambiguous_effect_scale" for f in data.findings))

    def test_ordinary_fold_change_accepted_with_explicit_confirmation(self):
        """Ambiguous fold change is accepted only when caller provides explicit confirmation override."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "linear_fc.tsv"
            p.write_text(
                "gene_id\tfold_change\tpvalue\tpadj\n"
                "ENSG1\t2.0\t0.01\t0.05\n",
                encoding="utf-8",
            )
            # Caller explicitly states this column represents log2fc
            data = read_deg_table(p, column_mapping={"log2fc": "fold_change"})
            self.assertTrue(data.is_valid)
            self.assertAlmostEqual(2.0, data.rows[0].log2_fold_change)


if __name__ == "__main__":
    unittest.main()
