"""Integration tests for route execution dispatch layer connecting router to quantification bridge."""

import gzip
from pathlib import Path
import tempfile
import unittest

from metadata.fixtures import canonical_bulk_study, reported, with_provenance
from src.ingestion.models import AssetStatus, RegisteredAsset, Route
from src.qc.validator import validate_bulk_counts
from src.quantification_bridge.adapter import convert_featurecounts_to_matrix
from src.quantification_bridge.manifest import APPROVED_PILOT_RUNS
from src.quantification_bridge.pairing import FastqPair, FastqPairingError, group_fastq_pairs
from src.routing import (
    RouteExecutionPlan,
    RouteExecutionResult,
    dispatch_route_execution,
    plan_route_execution,
    route_metadata_assets,
)
from tests.routing.fixtures import add_asset, routing_metadata, write_asset


FASTQ_R1 = "@read1/1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@read2/1\nTGCATGCATGCAAAAA\n+\nJJJJJJJJJJJJJJJJ\n"
FASTQ_R2 = "@read1/2\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@read2/2\nTTTTTGCATGCAAGCA\n+\nJJJJJJJJJJJJJJJJ\n"


def create_paired_fastq_study(directory: Path, sample_id: str = "sample_c1") -> dict:
    """Create canonical metadata and files for a paired FASTQ library."""
    r1_path = directory / f"{sample_id}_1.fastq.gz"
    r2_path = directory / f"{sample_id}_2.fastq.gz"

    with gzip.open(r1_path, "wt", encoding="utf-8") as f:
        f.write(FASTQ_R1)
    with gzip.open(r2_path, "wt", encoding="utf-8") as f:
        f.write(FASTQ_R2)

    doc = canonical_bulk_study()
    doc["sequencing_assays"][0]["read_layout"] = reported("paired_end")

    # Asset 1 (R1)
    a1 = doc["input_data_assets"][0]
    a1["asset_id"] = f"{sample_id}_r1_asset"
    a1["input_data_format"] = reported("FASTQ")
    a1["source_file_identifier_or_reference"] = reported(str(r1_path))
    a1.pop("gene_identifier_type", None)

    # Asset 2 (R2)
    a2 = dict(a1)
    a2["asset_id"] = f"{sample_id}_r2_asset"
    a2["source_file_identifier_or_reference"] = reported(str(r2_path))
    doc["input_data_assets"].append(a2)

    return with_provenance(doc)


class RouterDispatchIntegrationTests(unittest.TestCase):
    """Verify route execution planning, dispatch, pairing, and handoff contracts."""

    def test_fastq_route_planning_and_provenance(self):
        """Valid paired FASTQ input plans quantification_bridge execution with full provenance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata = create_paired_fastq_study(Path(tmpdir), sample_id="sample_c1")
            plan = plan_route_execution(metadata)

            self.assertEqual(Route.RAW_READ_PIPELINE, plan.route)
            self.assertEqual("READY_FOR_EXECUTION", plan.status)
            self.assertEqual("quantification_bridge", plan.executor_id)
            self.assertEqual("raw_counts", plan.target_format)
            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, plan.downstream_route)

            # Check provenance chain
            prov = plan.provenance_chain
            self.assertEqual("FASTQ", prov["detected_format"])
            self.assertEqual("RAW_READ_PIPELINE", prov["selected_route"])
            self.assertEqual("quantification_bridge", prov["raw_read_executor"])
            self.assertEqual("raw_counts", prov["emitted_format"])
            self.assertEqual("BULK_RAW_COUNT_PIPELINE", prov["downstream_route"])

    def test_raw_counts_route_planning(self):
        """Valid bulk integer matrix plans direct bulk_count_pipeline_runner execution."""
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_path = write_asset(Path(tmpdir) / "counts.tsv", "integer_matrix")
            metadata = routing_metadata(matrix_path, "raw_counts")
            plan = plan_route_execution(metadata)

            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, plan.route)
            self.assertEqual("READY_FOR_EXECUTION", plan.status)
            self.assertEqual("bulk_count_pipeline_runner", plan.executor_id)
            self.assertEqual("normalized_dataset", plan.target_format)
            self.assertIsNone(plan.downstream_route)

    def test_no_cross_route_contamination_fastq_cannot_enter_step6a(self):
        """Raw FASTQ cannot be parsed by Step 6A bulk raw-count validator."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq_path = write_asset(Path(tmpdir) / "reads_1.fastq.gz", "fastq_gz")
            doc = canonical_bulk_study()
            doc["input_data_assets"][0]["source_file_identifier_or_reference"] = reported(str(fastq_path))
            doc = with_provenance(doc)

            # Direct attempt to run Step 6A on FASTQ fails safely
            res = validate_bulk_counts(fastq_path, doc, "counts_asset", "bulk_assay")
            self.assertEqual("FAIL", res.status)

            # Router never maps FASTQ to BULK_RAW_COUNT_PIPELINE
            plan = plan_route_execution(doc)
            self.assertNotEqual(Route.BULK_RAW_COUNT_PIPELINE, plan.route)

    def test_no_cross_route_contamination_raw_counts_never_trigger_bridge(self):
        """Count matrix plan never selects the quantification bridge executor."""
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_path = write_asset(Path(tmpdir) / "counts.tsv", "integer_matrix")
            metadata = routing_metadata(matrix_path, "raw_counts")
            plan = plan_route_execution(metadata)

            self.assertNotEqual("quantification_bridge", plan.executor_id)
            self.assertEqual("bulk_count_pipeline_runner", plan.executor_id)

    def test_invalid_fastq_halts_before_execution(self):
        """Malformed FASTQ produces INVALID_ASSETS plan and dispatch halts."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bad_path = write_asset(Path(tmpdir) / "bad.fastq", "malformed_fastq")
            metadata = routing_metadata(bad_path, "FASTQ")
            plan = plan_route_execution(metadata)

            self.assertEqual("INVALID_ASSETS", plan.status)
            self.assertIsNone(plan.executor_id)

            # Dispatching halts with failure
            res = dispatch_route_execution(plan)
            self.assertFalse(res.is_success)
            self.assertEqual("INVALID_ASSETS", res.execution_state)

    def test_declared_detected_mismatch_preserves_manual_review(self):
        """Declared/detected mismatch routes to MANUAL_REVIEW_REQUIRED and halts dispatch."""
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_path = write_asset(Path(tmpdir) / "matrix.tsv", "processed_matrix")
            metadata = routing_metadata(matrix_path, "raw_counts")
            plan = plan_route_execution(metadata)

            self.assertEqual(Route.MANUAL_REVIEW, plan.route)
            self.assertEqual("MANUAL_REVIEW_REQUIRED", plan.status)
            self.assertIsNone(plan.executor_id)

            res = dispatch_route_execution(plan)
            self.assertFalse(res.is_success)
            self.assertEqual("MANUAL_REVIEW_REQUIRED", res.execution_state)

    def test_paired_read_grouping_deterministic(self):
        """group_fastq_pairs accurately pairs R1 and R2 regardless of input asset ordering."""
        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            p1 = write_asset(td / "s1_2.fastq.gz", "fastq_gz")
            p2 = write_asset(td / "s1_1.fastq.gz", "fastq_gz")
            p3 = write_asset(td / "s2_1.fastq.gz", "fastq_gz")
            p4 = write_asset(td / "s2_2.fastq.gz", "fastq_gz")

            # Feed in reverse / mixed order
            assets = (
                RegisteredAsset("a4", str(p4), "FASTQ", ("assay_s2",), ("bulk_rna_seq",)),
                RegisteredAsset("a1", str(p1), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
                RegisteredAsset("a3", str(p3), "FASTQ", ("assay_s2",), ("bulk_rna_seq",)),
                RegisteredAsset("a2", str(p2), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
            )

            pairs = group_fastq_pairs(assets)
            self.assertEqual(2, len(pairs))
            # Sort is deterministically s1, then s2
            self.assertEqual("s1", pairs[0].sample_id)
            self.assertEqual(p2, pairs[0].r1_path)
            self.assertEqual(p1, pairs[0].r2_path)

            self.assertEqual("s2", pairs[1].sample_id)
            self.assertEqual(p3, pairs[1].r1_path)
            self.assertEqual(p4, pairs[1].r2_path)

    def test_paired_read_integrity_orphan_mate_fails(self):
        """An orphan forward or reverse read halts grouping."""
        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            p1 = write_asset(td / "s1_1.fastq.gz", "fastq_gz")
            assets = (
                RegisteredAsset("a1", str(p1), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
            )
            with self.assertRaises(FastqPairingError) as ctx:
                group_fastq_pairs(assets)
            self.assertIn("Orphan mate detected", str(ctx.exception))

    def test_paired_read_integrity_duplicate_mate_fails(self):
        """Duplicate R1 reads for the same sample halt grouping."""
        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            p1 = write_asset(td / "s1_runA_1.fastq.gz", "fastq_gz")
            p2 = write_asset(td / "s1_runB_1.fastq.gz", "fastq_gz")
            p3 = write_asset(td / "s1_runA_2.fastq.gz", "fastq_gz")

            metadata = {
                "sequencing_assays": [
                    {"assay_id": "assay_s1", "sample_ids": ["s1"]}
                ]
            }

            assets = (
                RegisteredAsset("a1", str(p1), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
                RegisteredAsset("a2", str(p2), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
                RegisteredAsset("a3", str(p3), "FASTQ", ("assay_s1",), ("bulk_rna_seq",)),
            )
            with self.assertRaises(FastqPairingError) as ctx:
                group_fastq_pairs(assets, metadata=metadata)
            self.assertIn("Duplicate", str(ctx.exception))

    def test_sra_route_plans_and_dispatches_extraction(self):
        """SRA container input routes to RAW_READ_PIPELINE with fasterq-dump extraction adapter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sra_path = write_asset(Path(tmpdir) / "reads.sra", "sra")
            metadata = routing_metadata(sra_path, "SRA")
            plan = plan_route_execution(metadata)

            self.assertEqual(Route.RAW_READ_PIPELINE, plan.route)
            self.assertEqual("READY_FOR_EXECUTION", plan.status)
            self.assertEqual("quantification_bridge", plan.executor_id)
            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, plan.downstream_route)
            self.assertEqual("fasterq-dump", plan.provenance_chain.get("sra_extraction_adapter"))

    def test_synthetic_raw_read_output_satisfies_frozen_step6a_contract(self):
        """Synthetic featureCounts converted matrix satisfies the exact Step 6A bulk count contract."""
        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            # Create synthetic featureCounts native file
            fc_file = td / "synthetic_featureCounts.txt"
            fc_lines = [
                "# Program:featureCounts v2.0.6; -p --countReadPairs -B -C -s 2 -a gtf -o out bam1 bam2\n",
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/data/align/sample_c1_Aligned.sortedByCoord.out.bam\t/data/align/sample_c2_Aligned.sortedByCoord.out.bam\n",
                "ENSG00000000003\tchrX\t100\t500\t-\t401\t250\t310\n",
                "ENSG00000000005\tchrX\t1000\t2000\t+\t1001\t0\t5\n",
                "ENSG00000000419\tchr20\t500\t900\t-\t401\t1200\t1150\n",
            ]
            fc_file.write_text("".join(fc_lines), encoding="utf-8")

            # Convert using governed adapter
            out_matrix = td / "raw_counts.tsv"
            adapter_res = convert_featurecounts_to_matrix(
                raw_featurecounts_path=fc_file,
                output_matrix_path=out_matrix,
                approved_runs=("sample_c1", "sample_c2"),
            )
            self.assertTrue(out_matrix.exists())
            self.assertEqual(3, adapter_res["total_genes_emitted"])
            self.assertEqual(0, adapter_res["numerical_conservation_discrepancies"])

            # Verify that frozen Step 6A validator consumes this matrix cleanly
            doc = canonical_bulk_study()
            # Set up sample_c1 and sample_c2 in metadata
            doc["samples"] = [
                {
                    "sample_id": "sample_c1",
                    "experiment_id": "exp_001",
                    "condition_id": "cond_ctrl",
                    "biological_source_id": "source_001",
                    "organoid_context_id": "context_001",
                    "biological_replicate_id": reported("rep_1"),
                    "developmental_age_or_stage_at_collection": reported("Day 30"),
                    "metadata_validation_status": "valid",
                },
                {
                    "sample_id": "sample_c2",
                    "experiment_id": "exp_001",
                    "condition_id": "cond_ctrl",
                    "biological_source_id": "source_001",
                    "organoid_context_id": "context_001",
                    "biological_replicate_id": reported("rep_2"),
                    "developmental_age_or_stage_at_collection": reported("Day 30"),
                    "metadata_validation_status": "valid",
                },
            ]
            doc["sequencing_assays"][0]["sample_ids"] = ["sample_c1", "sample_c2"]
            doc["input_data_assets"][0]["source_file_identifier_or_reference"] = reported(str(out_matrix))
            doc = with_provenance(doc)

            qc_res = validate_bulk_counts(out_matrix, doc, "counts_asset", "bulk_assay")
            # Step 6A should successfully process the matrix rows and samples
            self.assertEqual(3, qc_res.number_of_genes)
            self.assertEqual(2, qc_res.number_of_samples)
            self.assertIn(qc_res.status, {"PASS", "PASS_WITH_WARNINGS"})

    def test_pilot_12_run_manifest_invariance(self):
        """PRJNA802072 pilot manifest remains strictly 12 non-contiguous runs."""
        self.assertEqual(12, len(APPROVED_PILOT_RUNS))
        controls = [r for r in APPROVED_PILOT_RUNS if r.condition_group == "control"]
        treatments = [r for r in APPROVED_PILOT_RUNS if r.condition_group == "treatment"]
        self.assertEqual(6, len(controls))
        self.assertEqual(6, len(treatments))
        # Ensure non-contiguous runs
        run_ids = [r.run_accession for r in APPROVED_PILOT_RUNS]
        self.assertEqual(12, len(set(run_ids)))
        self.assertIn("SRR17834920", run_ids)
        self.assertIn("SRR17834900", run_ids)


if __name__ == "__main__":
    unittest.main()
