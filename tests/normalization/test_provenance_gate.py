"""Integration and unit tests for Phase U3.2 Processing Provenance Gate."""

from pathlib import Path
import tempfile
import unittest

from tests.metadata.fixtures import canonical_bulk_study, reported, with_provenance
from src.gene_harmonization.models import GeneMappingStatus
from src.ingestion.processed_matrix_reader import ProcessedExpressionMatrix
from src.normalization.provenance_gate import (
    ExpressionUnit,
    NormalizationMethod,
    ProcessedExpressionDataset,
    ProcessedNormalizationMethod,
    ProcessingProvenance,
    ProvenanceGateResult,
    TransformationScale,
    evaluate_provenance_gate,
    ingest_processed_expression_matrix,
    parse_processing_provenance,
)
from src.qc.result import Severity, Status
from src.qc.validator import validate_bulk_counts


def _make_metadata(
    samples: list[dict],
    processing_provenance: dict | None = None,
    study_id: str = "study_u3_test",
) -> dict:
    doc = canonical_bulk_study()
    doc["studies"][0]["study_id"] = study_id
    doc["samples"] = samples
    sample_ids = [s["sample_id"] for s in samples]
    doc["sequencing_assays"][0]["sample_ids"] = sample_ids
    doc["input_data_assets"][0]["assay_ids"] = [doc["sequencing_assays"][0]["assay_id"]]
    if processing_provenance is not None:
        doc["processing_provenance"] = processing_provenance
        doc["input_data_assets"][0]["processing_provenance"] = processing_provenance
    return with_provenance(doc)


class ProcessingProvenanceGateTests(unittest.TestCase):
    """Test Processing Provenance Gate: units, scales, no-guessing rules, compatibility, and DESeq2 prevention."""

    def test_valid_tpm_matrix_with_documented_provenance(self):
        """1. Valid TPM matrix with linear scale and documented software pipeline qualifies for limma-trend."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("ENSG00000141510", "ENSG00000012048"),
            sample_ids=("s1", "s2"),
            columns=((15.2, 4.8), (18.1, 5.2)),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "s1", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "s2", "condition_id": "c1", "biological_replicate_id": reported("rep2")},
            ],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
                "software_pipeline": "Salmon v1.10.0",
                "reference_annotation": "GENCODE v44",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)

        self.assertEqual(Status.PASS, res.status)
        self.assertTrue(res.is_inferentially_eligible)
        self.assertEqual(ExpressionUnit.TPM, res.provenance.expression_unit)
        self.assertEqual(TransformationScale.LINEAR, res.provenance.transformation_scale)
        self.assertEqual("limma_trend", res.target_analysis_engine)
        self.assertFalse(res.deseq2_eligible)
        self.assertIn("DESeq2", res.blocked_analysis_engines)

    def test_valid_fpkm_matrix_with_documented_provenance(self):
        """2. Valid FPKM matrix with linear scale qualifies for downstream analysis."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("ENSG00000141510", "ENSG00000012048"),
            sample_ids=("s1", "s2"),
            columns=((25.4, 8.1), (28.9, 9.0)),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "s1", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "s2", "condition_id": "c1", "biological_replicate_id": reported("rep2")},
            ],
            processing_provenance={
                "expression_unit": "FPKM",
                "transformation_scale": "linear",
                "normalization_method": "fpkm",
                "software_pipeline": "StringTie v2.2.1",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)

        self.assertEqual(Status.PASS, res.status)
        self.assertTrue(res.is_inferentially_eligible)
        self.assertEqual(ExpressionUnit.FPKM, res.provenance.expression_unit)
        self.assertFalse(res.provenance.is_log_transformed)

    def test_valid_cpm_and_logcpm_matrix(self):
        """3. Valid CPM and logCPM matrices with known transformation scale."""
        # Case A: CPM linear
        matrix_cpm = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((120.5, 45.2),),
        )
        meta_cpm = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "CPM",
                "transformation_scale": "linear",
                "normalization_method": "cpm",
                "software_pipeline": "edgeR::cpm",
            },
        )
        res_cpm = evaluate_provenance_gate(matrix_cpm, meta_cpm)
        self.assertEqual(Status.PASS, res_cpm.status)
        self.assertFalse(res_cpm.provenance.is_log_transformed)

        # Case B: logCPM (inherently log2 scale)
        matrix_logcpm = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((6.91, 5.50),),
        )
        meta_logcpm = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "logCPM",
                "normalization_method": "cpm",
                "software_pipeline": "edgeR::cpm(log=TRUE)",
            },
        )
        res_logcpm = evaluate_provenance_gate(matrix_logcpm, meta_logcpm)
        self.assertEqual(Status.PASS, res_logcpm.status)
        self.assertEqual(ExpressionUnit.LOG_CPM, res_logcpm.provenance.expression_unit)
        self.assertEqual(TransformationScale.LOG2, res_logcpm.provenance.transformation_scale)
        self.assertTrue(res_logcpm.provenance.is_log_transformed)

    def test_valid_log_transformed_tpm_with_pseudocount(self):
        """Log2-transformed TPM (log2(TPM+1)) matrix with documented offset."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((4.15, 0.0),),  # 0.0 corresponds to log2(0 + 1)
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "log2(x+1)",
                "pseudocount": 1.0,
                "normalization_method": "tpm",
                "software_pipeline": "RSEM v1.3.3",
            },
        )
        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.PASS, res.status)
        self.assertTrue(res.provenance.is_log_transformed)
        self.assertEqual(1.0, res.provenance.pseudocount)

    def test_valid_size_factor_normalized_counts(self):
        """4. Normalized count matrix with documented size factor normalization method."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1", "s2"),
            columns=((150.25, 20.10), (145.80, 22.40)),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "s1", "condition_id": "c1"},
                {"sample_id": "s2", "condition_id": "c1"},
            ],
            processing_provenance={
                "expression_unit": "normalized_counts",
                "transformation_scale": "linear",
                "normalization_method": "deseq2_size_factors",
                "software_pipeline": "DESeq2::counts(normalized=TRUE)",
            },
        )
        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.PASS, res.status)
        self.assertEqual(ExpressionUnit.NORMALIZED_COUNTS, res.provenance.expression_unit)
        self.assertEqual(ProcessedNormalizationMethod.DESEQ2_SIZE_FACTORS, res.provenance.normalization_method)

    def test_unknown_expression_unit_blocks_inferential_processing(self):
        """5A. Undocumented expression unit blocks processing with INSUFFICIENT_PROCESSING_PROVENANCE."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("s1",),
            columns=((12.5,),),
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                # Missing expression_unit
                "transformation_scale": "linear",
                "normalization_method": "none",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.FAIL, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertIsNone(res.target_analysis_engine)
        self.assertTrue(
            any(
                f.rule_id == "INSUFFICIENT_PROCESSING_PROVENANCE"
                and f.path == "expression_unit"
                and f.severity == Severity.ERROR
                for f in res.findings
            )
        )

    def test_unknown_transformation_scale_blocks_no_guessing(self):
        """5B. Undocumented transformation scale is NOT guessed; blocks with INSUFFICIENT_PROCESSING_PROVENANCE."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("s1",),
            columns=((12.5,),),  # Could be linear or log2; system MUST NOT GUESS
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "TPM",
                # transformation_scale is omitted!
                "normalization_method": "tpm",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.FAIL, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertTrue(
            any(
                f.rule_id == "INSUFFICIENT_PROCESSING_PROVENANCE"
                and f.path == "transformation_scale"
                and f.severity == Severity.ERROR
                for f in res.findings
            )
        )

    def test_negative_values_in_linear_scale_rejected(self):
        """8A. Negative values in declared linear expression scale emit NEGATIVE_VALUES_IN_LINEAR_SCALE ERROR."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((12.5, -0.45),),  # Negative value in linear scale!
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.FAIL, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertTrue(any(f.rule_id == "NEGATIVE_VALUES_IN_LINEAR_SCALE" and f.severity == Severity.ERROR for f in res.findings))

    def test_negative_values_in_offset_log_scale_rejected(self):
        """8B. Negative values in log2(x + 1) scale emit NEGATIVE_VALUES_IN_OFFSET_LOG_SCALE ERROR."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((3.5, -0.1),),  # Negative value in log2(x+1) for x>=0 is impossible
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "log2(x+1)",
                "pseudocount": 1.0,
                "normalization_method": "tpm",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.FAIL, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertTrue(any(f.rule_id == "NEGATIVE_VALUES_IN_OFFSET_LOG_SCALE" and f.severity == Severity.ERROR for f in res.findings))

    def test_non_finite_values_rejected(self):
        """8C. Non-finite values (NaN / Inf) emit NON_FINITE_EXPRESSION_VALUES ERROR."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1",),
            columns=((10.0, float("nan")),),
        )
        meta = _make_metadata(
            samples=[{"sample_id": "s1", "condition_id": "c1"}],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
            },
        )

        res = evaluate_provenance_gate(matrix, meta)
        self.assertEqual(Status.FAIL, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertTrue(any(f.rule_id == "NON_FINITE_EXPRESSION_VALUES" for f in res.findings))

    def test_prevention_of_accidental_routing_into_deseq2(self):
        """9. Processed continuous expression data can NEVER pass into DESeq2 raw-count validator."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\ts1\ts2\n")
            f.write("G1\t12.45\t0.50\n")
            f.write("G2\t3.10\t18.90\n")
            tmp_path = Path(f.name)

        meta = _make_metadata(
            samples=[
                {"sample_id": "s1", "condition_id": "c1"},
                {"sample_id": "s2", "condition_id": "c1"},
            ]
        )

        try:
            # Feeding continuous processed data into raw-count Step 6A validator
            raw_qc_res = validate_bulk_counts(
                tmp_path,
                meta,
                asset_id="asset_test",
                assay_id=meta["sequencing_assays"][0]["assay_id"],
            )
            # MUST FAIL because values are not integers!
            self.assertEqual(Status.FAIL, raw_qc_res.status)
            self.assertTrue(any(f.rule_id == "integer_like_counts_required" for f in raw_qc_res.findings))
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_end_to_end_processed_expression_intake_and_harmonization(self):
        """10. End-to-end intake, QC, provenance gating, and gene reference harmonization."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        self.assertTrue(fixture_path.exists())

        meta = _make_metadata(
            samples=[
                {"sample_id": "GSM7321024", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321026", "condition_id": "cond_vpa200", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321027", "condition_id": "cond_vpa400", "biological_replicate_id": reported("rep1")},
            ],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
                "software_pipeline": "Kallisto v0.48.0",
                "reference_annotation": "Ensembl Release 112",
            },
            study_id="GSE232218",
        )
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]

        dataset = ingest_processed_expression_matrix(fixture_path, meta)

        self.assertEqual(Status.PASS, dataset.status)
        self.assertTrue(dataset.is_valid_intake)
        self.assertTrue(dataset.is_provenance_valid)
        self.assertTrue(dataset.is_u3_3_evaluable)
        # Disentangled eligibility: GSE232218 has n=1 per condition, so statistical differential inference is blocked
        self.assertFalse(dataset.is_statistically_eligible)
        self.assertFalse(dataset.is_inferentially_eligible)
        self.assertIsNotNone(dataset.statistical_eligibility_reason)
        self.assertIn("fewer than 2 biological replicates", dataset.statistical_eligibility_reason.lower())
        self.assertIsNone(dataset.target_analysis_engine)
        self.assertIn("DESeq2", dataset.blocked_analysis_engines)

        # Matrix dimensions
        self.assertEqual(15, dataset.matrix.number_of_genes)
        self.assertEqual(3, dataset.matrix.number_of_samples)

        # Exact numerical preservation: GSM7321024 TP53 == 42.1582
        self.assertEqual(42.1582, dataset.matrix.get_value("ENSG00000141510", "GSM7321024"))

        # Gene Harmonization against Ensembl 112
        self.assertEqual(15, len(dataset.harmonized_genes))
        tp53_gene = next(g for g in dataset.harmonized_genes if g.original_gene_id == "ENSG00000141510")
        self.assertEqual("ENSG00000141510", tp53_gene.canonical_gene_id)
        self.assertEqual("TP53", tp53_gene.approved_symbol)

        # Law of Conservation of the Gene Universe
        self.assertTrue(dataset.is_gene_universe_conserved)
        self.assertEqual(15, dataset.total_input_genes_count)
        self.assertEqual(
            dataset.total_input_genes_count,
            dataset.eligible_canonical_genes_count + dataset.excluded_genes_count,
        )

    def test_replicated_conditions_are_statistically_and_inferentially_eligible(self):
        """11. Replicated experimental conditions (n >= 2) qualify for limma-trend differential inference."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = _make_metadata(
            samples=[
                {"sample_id": "GSM7321024", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321026", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep2")},
                {"sample_id": "GSM7321027", "condition_id": "cond_treated", "biological_replicate_id": reported("rep1")},
            ],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
                "software_pipeline": "Kallisto v0.48.0",
                "reference_annotation": "Ensembl Release 112",
            },
            study_id="GSE232218",
        )
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]

        # First case: cond_treated has only 1 replicate -> ineligible
        dataset1 = ingest_processed_expression_matrix(fixture_path, meta)
        self.assertFalse(dataset1.is_statistically_eligible)
        self.assertFalse(dataset1.is_inferentially_eligible)

        # Second case: assign GSM7321027 to cond_vehicle as rep3 -> condition has 3 replicates
        meta["samples"][2]["condition_id"] = "cond_vehicle"
        meta["samples"][2]["biological_replicate_id"] = reported("rep3")
        dataset2 = ingest_processed_expression_matrix(fixture_path, meta)
        self.assertTrue(dataset2.is_statistically_eligible)
        self.assertTrue(dataset2.is_inferentially_eligible)
        self.assertEqual("limma_trend", dataset2.target_analysis_engine)

    def test_production_gene_reference_unavailable_preserves_data_with_review_finding(self):
        """12. Missing production gene reference does NOT silently fallback; preserves data with NEEDS_REVIEW finding."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = _make_metadata(
            samples=[
                {"sample_id": "GSM7321024", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321026", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep2")},
                {"sample_id": "GSM7321027", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep3")},
            ],
            processing_provenance={
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
            },
            study_id="GSE232218",
        )
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]

        # Provide a config pointing to a nonexistent reference file
        bad_config = {
            "config_version": "0.1.0",
            "stage_name": "gene_harmonization",
            "scope": {"species": "Homo sapiens", "taxonomy_id": 9606, "modality": "bulk_rna_seq", "input_format": "raw_counts"},
            "reference": {
                "reference_id": "missing_ref",
                "reference_version": "1",
                "genome_build": "GRCh38",
                "canonical_source_authority": "Ensembl",
                "file_path": "nonexistent/path/to/reference.tsv",
            },
            "supported_identifier_types": {},
            "structural_rules": {},
            "coverage_thresholds": {},
        }

        dataset = ingest_processed_expression_matrix(fixture_path, meta, config=bad_config)

        # Must report NEEDS_REVIEW due to reference unavailability
        self.assertEqual(Status.NEEDS_REVIEW, dataset.status)
        self.assertFalse(dataset.is_inferentially_eligible)
        self.assertTrue(any(f.rule_id == "gene_reference_unavailable" and f.severity == Severity.REVIEW for f in dataset.gate_result.findings))
        # Original gene IDs must still be preserved unmapped
        self.assertEqual(15, len(dataset.harmonized_genes))
        self.assertTrue(all(g.mapping_status == GeneMappingStatus.UNMAPPED for g in dataset.harmonized_genes))


if __name__ == "__main__":
    unittest.main()
