"""Unit, integration, and real-data validation tests for published differential-expression adapter (Phase U3.1)."""

import math
from pathlib import Path
import tempfile
import unittest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.differential_expression.published_adapter import (
    derive_published_response_direction,
    ingest_published_deg_table,
)
from src.features.cross_study.builder import build_records_for_contrast
from src.features.cross_study.models import (
    EstimateStatus,
    ResponseDirection,
)
from src.response_builder.models import ContrastStatus


class PublishedDegAdapterTests(unittest.TestCase):
    """Verify published DEG table adaptation, reference harmonization, and Stage 8A.1 integration."""

    def test_descriptive_response_direction_semantics(self):
        """Verify strict Stage 8A.1 descriptive response direction semantics:
        UP (log2FC > 0), DOWN (log2FC < 0), NO_DIRECTION (log2FC == 0), UNDETERMINED (unestimable).
        """
        self.assertEqual(ResponseDirection.UP, derive_published_response_direction(1.5))
        self.assertEqual(ResponseDirection.UP, derive_published_response_direction(0.0001))
        self.assertEqual(ResponseDirection.DOWN, derive_published_response_direction(-0.8))
        self.assertEqual(ResponseDirection.DOWN, derive_published_response_direction(-5.2))
        self.assertEqual(ResponseDirection.NO_DIRECTION, derive_published_response_direction(0.0))
        self.assertEqual(ResponseDirection.UNDETERMINED, derive_published_response_direction(None))
        self.assertEqual(ResponseDirection.UNDETERMINED, derive_published_response_direction(float("nan")))
        self.assertEqual(ResponseDirection.UNDETERMINED, derive_published_response_direction(float("inf")))

    def test_ingest_synthetic_deseq2_table(self):
        """Synthetic DESeq2 table with Ensembl IDs maps to canonical reference and preserves values."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "deg_synthetic.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\tgene_name\n"
                "ENSG00000101017\t3.397\t2.62e-106\t7.21e-102\tCD40\n"
                "ENSG00000130876\t-1.850\t0.00045\t0.0032\tSLC7A10\n"
                "ENSG00000000000\t0.500\t0.01\t0.05\tFAKE_GENE\n",
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "agent_name": "valproic acid",
                "concentration_or_dose": 1.0,
                "concentration_or_dose_unit": "millimolar",
                "exposure_duration": 29.0,
                "exposure_duration_unit": "days",
                "treatment_replicate_count": 8,
                "control_replicate_count": 8,
            }

            contrast_res, table_data = ingest_published_deg_table(p, meta)

            self.assertEqual(ContrastStatus.ELIGIBLE, contrast_res.status)
            self.assertEqual(3, contrast_res.total_input_genes_count)
            self.assertEqual(2, contrast_res.eligible_canonical_genes_count)
            self.assertEqual(1, contrast_res.excluded_genes_count)
            self.assertTrue(contrast_res.is_gene_universe_conserved)

            # Check unmapped gene was audited into excluded_genes, not silently dropped
            self.assertEqual("ENSG00000000000", contrast_res.excluded_genes[0].original_gene_id)
            self.assertEqual("UNMAPPED", contrast_res.excluded_genes[0].mapping_status)

            # Check canonical results
            cd40 = contrast_res.gene_results[0]
            self.assertEqual("ENSG00000101017", cd40.canonical_gene_id)
            self.assertEqual("CD40", cd40.approved_symbol)
            self.assertAlmostEqual(3.397, cd40.log2_fold_change)
            self.assertAlmostEqual(2.62e-106, cd40.p_value)
            self.assertAlmostEqual(7.21e-102, cd40.adjusted_p_value_bh)

            # Check zero data fabrication: missing standard error remains None
            self.assertIsNone(cd40.lfc_standard_error)
            self.assertIsNone(cd40.wald_statistic)

    def test_direction_separated_from_statistical_significance(self):
        """CRITICAL: Non-significant genes (padj > 0.05) with positive log2FC have direction UP, NOT UNDETERMINED."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "nonsig.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "ENSG00000101017\t0.15\t0.45\t0.78\n",  # Non-significant, but positive log2FC
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "agent_name": "lead",
                "treatment_replicate_count": 3,
                "control_replicate_count": 3,
            }

            contrast_res, _ = ingest_published_deg_table(p, meta)
            records = build_records_for_contrast(contrast_res)
            self.assertEqual(1, len(records))

            rec = records[0]
            # Direction is strictly descriptive based on sign of log2FC
            self.assertEqual(ResponseDirection.UP, rec.response_direction)
            self.assertAlmostEqual(0.15, rec.log2_fold_change)
            self.assertAlmostEqual(0.78, rec.adjusted_p_value)

    def test_unexplained_missing_pvalues_remain_unexplained_not_cooks_outlier(self):
        """Missing p-values without source-verified outlier reason remain unexplained missing statistics."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "missing_p.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "ENSG00000101017\t4.5\tNA\tNA\n",
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "agent_name": "rotenone",
                "treatment_replicate_count": 3,
                "control_replicate_count": 3,
            }

            contrast_res, _ = ingest_published_deg_table(p, meta)
            self.assertEqual(1, len(contrast_res.gene_results))
            g = contrast_res.gene_results[0]
            # Must NOT be labeled OUTLIER when no source reason exists
            self.assertEqual("MISSING_STATISTICS", g.result_status)
            self.assertIsNone(g.p_value)

            records = build_records_for_contrast(contrast_res)
            self.assertEqual(1, len(records))
            # Frozen Stage 8A.1 model classifies missing estimates without outlier label as INVALID_STATISTIC
            self.assertEqual(EstimateStatus.INVALID_STATISTIC, records[0].estimate_status)
            self.assertEqual(ResponseDirection.UP, records[0].response_direction)

    def test_source_established_outlier_preserved_as_cooks_outlier(self):
        """When source table explicitly establishes an outlier reason, it is mapped to COOKS_OUTLIER."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "source_outlier.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\tstatus\n"
                "ENSG00000101017\t4.5\tNA\tNA\tcooks_outlier\n",
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "agent_name": "rotenone",
                "treatment_replicate_count": 3,
                "control_replicate_count": 3,
            }

            contrast_res, _ = ingest_published_deg_table(p, meta)
            self.assertEqual(1, len(contrast_res.gene_results))
            self.assertEqual("OUTLIER", contrast_res.gene_results[0].result_status)

            records = build_records_for_contrast(contrast_res)
            self.assertEqual(1, len(records))
            self.assertEqual(EstimateStatus.COOKS_OUTLIER, records[0].estimate_status)

    def test_duplicate_canonical_collisions_audited(self):
        """Duplicate rows mapping to the same canonical Ensembl ID are audited as collisions without silent averaging."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "collisions.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "ENSG00000101017\t3.0\t0.001\t0.01\n"
                "ENSG00000101017\t3.5\t0.002\t0.02\n",
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "treatment_replicate_count": 3,
                "control_replicate_count": 3,
            }

            contrast_res, _ = ingest_published_deg_table(p, meta)
            self.assertEqual(1, contrast_res.eligible_canonical_genes_count)
            self.assertEqual(1, contrast_res.excluded_genes_count)
            self.assertTrue(contrast_res.is_gene_universe_conserved)

            collision = contrast_res.excluded_genes[0]
            self.assertEqual("COLLISION", collision.mapping_status)
            self.assertEqual("ENSG00000101017", collision.canonical_gene_id)

    def test_under_replicated_contrast_flagged(self):
        """Contrast with only 1 biological replicate is flagged for review rather than ELIGIBLE."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "single_rep.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\n"
                "ENSG00000101017\t2.0\t0.01\t0.05\n",
                encoding="utf-8",
            )
            meta = {
                "dataset": {"study_id": "test_study"},
                "treatment_replicate_count": 1,
                "control_replicate_count": 1,
            }

            contrast_res, _ = ingest_published_deg_table(p, meta)
            self.assertEqual(ContrastStatus.NEEDS_REVIEW, contrast_res.status)
            self.assertFalse(contrast_res.is_inferentially_eligible)

    def test_real_published_dataset_vpa_dorsal_forebrain_organoid_tsv(self):
        """Real published DEG validation: Human dorsal forebrain organoids exposed to valproate (PMC13441967).
        Reconciled experimental metadata: 1 mM VPA from Day 21 to Day 50 (29 days) across 4 differentiation batches (N=8).
        """
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        self.assertTrue(fixture_path.is_file(), f"Real fixture missing at {fixture_path}")

        meta = {
            "dataset": {
                "study_id": "E-MTAB-15191",
                "cohort_id": "cohort_dorsal_forebrain_organoids_bioni010c",
            },
            "agent_name": "valproic acid",
            "agent_identifier": "CHEBI:39867",
            "concentration_or_dose": 1.0,
            "concentration_or_dose_unit": "millimolar",
            "exposure_start_time_or_stage": "Day 21 of differentiation",
            "developmental_age_or_stage_at_exposure": "Day 21 of differentiation",
            "treatment_collection_age_or_stage": "Day 50",
            "exposure_duration": 29.0,
            "exposure_duration_unit": "days",
            "treatment_replicate_count": 8,
            "control_replicate_count": 8,
            "treatment_sample_ids": ["AB66V1", "AB66V2", "AB67V1", "AB67V2", "AB69V1", "AB69V2", "AB70V1", "AB70V2"],
            "control_sample_ids": ["AB66C1", "AB66C2", "AB67C1", "AB67C2", "AB69C1", "AB69C2", "AB70C1", "AB70C2"],
            "treatment_condition_id": "VPA_1mM_D21_D50",
            "control_condition_id": "Vehicle_Control_D50",
            "published_analysis_method": "Novogene (unspecified differential expression pipeline; Fig 2B reports t-test padj, Table S2 provides log2FoldChange, pvalue, padj)",
        }

        contrast_res, table_data = ingest_published_deg_table(fixture_path, meta)

        self.assertEqual(ContrastStatus.ELIGIBLE, contrast_res.status)
        self.assertTrue(contrast_res.is_inferentially_eligible)
        self.assertEqual(500, contrast_res.total_input_genes_count)
        self.assertTrue(contrast_res.eligible_canonical_genes_count >= 490)
        self.assertTrue(contrast_res.is_gene_universe_conserved)

        # Check provenance
        prov = contrast_res.r_environment_info
        self.assertEqual("published_deg_table_ingestion_v2_u3", prov["source_pipeline"])
        self.assertEqual("Novogene (unspecified differential expression pipeline; Fig 2B reports t-test padj, Table S2 provides log2FoldChange, pvalue, padj)", prov["analysis_method"])
        self.assertEqual(500, prov["total_rows_parsed"])

        # Check known published marker genes from the paper
        results_by_sym = {g.approved_symbol: g for g in contrast_res.gene_results if g.approved_symbol}

        # Top upregulated gene: CD40
        self.assertIn("CD40", results_by_sym)
        cd40 = results_by_sym["CD40"]
        self.assertEqual("ENSG00000101017", cd40.canonical_gene_id)
        self.assertAlmostEqual(3.39707825, cd40.log2_fold_change, places=4)
        self.assertAlmostEqual(2.62e-106, cd40.p_value, places=100)
        self.assertAlmostEqual(7.21e-102, cd40.adjusted_p_value_bh, places=95)

        # Other prominent neural progenitor / developmental markers from study
        for marker in ("SLC7A10", "ATP8B3", "DLX4", "BMP8B"):
            self.assertIn(marker, results_by_sym)
            m = results_by_sym[marker]
            self.assertTrue(m.log2_fold_change > 0)
            self.assertTrue(m.adjusted_p_value_bh < 0.05)

        # Seamless conversion to Stage 8A.1 ConditionGeneResponseRegistry records
        records = build_records_for_contrast(contrast_res)
        self.assertEqual(contrast_res.eligible_canonical_genes_count, len(records))

        # Check first record in Stage 8A.1 feature space
        r0 = records[0]
        self.assertEqual("ENSG00000101017", r0.canonical_ensembl_gene_id)
        self.assertEqual("CD40", r0.gene_symbol)
        self.assertEqual("valproic acid", r0.compound_name)
        self.assertEqual(1.0, r0.dose_value)
        self.assertEqual("millimolar", r0.dose_unit)
        self.assertEqual(29.0, r0.exposure_duration)
        self.assertEqual("days", r0.exposure_duration_unit)
        self.assertEqual(ResponseDirection.UP, r0.response_direction)
        self.assertIsNone(r0.log2_fold_change_se)  # Not invented!

    def test_real_published_dataset_vpa_dorsal_forebrain_organoid_xlsx(self):
        """Real published DEG validation via Excel (.xlsx) workbook with reconciled metadata."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.xlsx")
        self.assertTrue(fixture_path.is_file(), f"Real XLSX fixture missing at {fixture_path}")

        meta = {
            "dataset": {"study_id": "E-MTAB-15191"},
            "agent_name": "valproic acid",
            "concentration_or_dose": 1.0,
            "concentration_or_dose_unit": "millimolar",
            "exposure_duration": 29.0,
            "exposure_duration_unit": "days",
            "treatment_replicate_count": 8,
            "control_replicate_count": 8,
        }

        contrast_res, table_data = ingest_published_deg_table(
            fixture_path,
            meta,
            sheet_name="VPAvsCTRL_all_genes",
        )

        self.assertEqual(ContrastStatus.ELIGIBLE, contrast_res.status)
        self.assertEqual(500, contrast_res.total_input_genes_count)
        self.assertEqual("VPAvsCTRL_all_genes", table_data.sheet_name)
        self.assertTrue(contrast_res.is_gene_universe_conserved)

        # Check CD40 values in Excel format
        results_by_sym = {g.approved_symbol: g for g in contrast_res.gene_results if g.approved_symbol}
        self.assertIn("CD40", results_by_sym)
        self.assertAlmostEqual(3.39707825, results_by_sym["CD40"].log2_fold_change, places=4)

    def test_full_published_deg_table_38870_genes(self):
        """Test full published DEG table (Table S2, 38,870 genes) ensuring zero silent drop and complete conservation."""
        fixture_path = Path("tests/fixtures/Table_S2.xlsx")
        if not fixture_path.is_file():
            self.skipTest("tests/fixtures/Table_S2.xlsx not present")

        meta = {
            "dataset": {
                "study_id": "E-MTAB-15191",
                "cohort_id": "cohort_dorsal_forebrain_organoids_bioni010c",
            },
            "agent_name": "valproic acid",
            "agent_identifier": "CHEBI:39867",
            "concentration_or_dose": 1.0,
            "concentration_or_dose_unit": "millimolar",
            "exposure_duration": 29.0,
            "exposure_duration_unit": "days",
            "treatment_replicate_count": 8,
            "control_replicate_count": 8,
            "treatment_sample_ids": ["AB66V1", "AB66V2", "AB67V1", "AB67V2", "AB69V1", "AB69V2", "AB70V1", "AB70V2"],
            "control_sample_ids": ["AB66C1", "AB66C2", "AB67C1", "AB67C2", "AB69C1", "AB69C2", "AB70C1", "AB70C2"],
            "treatment_condition_id": "VPA_1mM_D21_to_D50",
            "control_condition_id": "Vehicle_Control_D50",
            "published_analysis_method": "Novogene (unspecified differential expression pipeline; Fig 2B reports t-test padj, Table S2 provides log2FoldChange, pvalue, padj)",
        }

        contrast_res, table_data = ingest_published_deg_table(
            fixture_path,
            meta,
            sheet_name="VPAvsCTRL_all_genes",
        )

        self.assertEqual(ContrastStatus.ELIGIBLE, contrast_res.status)
        self.assertEqual(38870, contrast_res.total_input_genes_count)
        self.assertEqual(38455, contrast_res.eligible_canonical_genes_count)
        self.assertEqual(415, contrast_res.excluded_genes_count)
        self.assertTrue(contrast_res.is_gene_universe_conserved)

        # Check conversion to Stage 8A.1 feature records for all 38,455 canonical genes
        records = build_records_for_contrast(contrast_res)
        self.assertEqual(38455, len(records))


    def test_default_analysis_method_is_unspecified_not_deseq2(self):
        """When published_analysis_method is not supplied, it defaults to 'unspecified', never 'DESeq2'."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = {
            "dataset": {"study_id": "E-MTAB-15191"},
            "agent_name": "valproic acid",
        }
        contrast_res, _ = ingest_published_deg_table(fixture_path, meta)
        prov = contrast_res.r_environment_info
        self.assertEqual("unspecified", prov["analysis_method"])
        self.assertNotIn("DESeq2", prov["analysis_method"])
        self.assertEqual("unspecified", contrast_res.design_formula)

if __name__ == "__main__":
    unittest.main()
