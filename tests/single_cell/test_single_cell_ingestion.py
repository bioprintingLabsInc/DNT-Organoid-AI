"""Acceptance and integration tests for Phase U4: Single-Cell / Single-Nucleus RNA-seq Ingestion & Pseudobulk."""

import gzip
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse

from src.ingestion.models import Route
from src.qc.matrix import CountMatrix, read_matrix
from src.qc.result import Status as QCStatus
from src.qc.validator import validate_bulk_counts
from src.routing import dispatch_route_execution, plan_route_execution
from src.single_cell import (
    CellMetadata,
    PseudobulkConfig,
    PseudobulkDataset,
    SingleCellMatrix,
    SingleCellQCConfig,
    SingleCellQCResult,
    UnsupportedSingleCellFormatException,
    aggregate_pseudobulk,
    detect_and_reject_rds,
    evaluate_single_cell_qc,
    read_10x_matrix_market,
    read_h5ad,
    run_single_cell_ingestion_pipeline,
)
from tests.metadata.fixtures import canonical_bulk_study, reported, with_provenance
from tests.routing.fixtures import routing_metadata


class PhaseU4SingleCellIngestionTests(unittest.TestCase):
    """Verify Phase U4 single-cell / single-nucleus ingestion, QC, pseudobulk, and downstream handoff."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def _create_synthetic_h5ad(
        self,
        n_cells: int = 50,
        n_genes: int = 20,
        is_sparse: bool = True,
    ) -> Path:
        """Create a valid synthetic AnnData H5AD file."""
        h5_path = self.dir_path / "synthetic_test.h5ad"

        gene_ids = [f"ENSG{100000 + i:011d}" for i in range(n_genes)]
        gene_symbols = [f"GENE_{i}" for i in range(n_genes)]
        gene_symbols[0] = "MT-CO1"
        gene_symbols[1] = "MT-ND1"

        cell_ids = [f"CELL_{j:04d}-1" for j in range(n_cells)]
        donors = [f"donor_{j % 3 + 1}" for j in range(n_cells)]  # 3 donors
        cell_types = ["cortical_neuron" if j % 2 == 0 else "radial_glia" for j in range(n_cells)]
        conditions = ["control" if j % 4 < 2 else "vpa_treated" for j in range(n_cells)]

        # Synthetic counts: shape (n_cells, n_genes) in standard AnnData
        rng = np.random.default_rng(42)
        raw_counts = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.int64)

        with h5py.File(h5_path, "w") as f:
            if is_sparse:
                csr = sparse.csr_matrix(raw_counts)
                x_group = f.create_group("X")
                x_group.attrs["encoding-type"] = "csr_matrix"
                x_group.attrs["shape"] = csr.shape
                x_group.create_dataset("data", data=csr.data)
                x_group.create_dataset("indices", data=csr.indices)
                x_group.create_dataset("indptr", data=csr.indptr)
            else:
                f.create_dataset("X", data=raw_counts)

            # obs group
            obs_group = f.create_group("obs")
            obs_group.attrs["_index"] = "cell_id"
            obs_group.create_dataset("cell_id", data=np.array(cell_ids, dtype=object))
            obs_group.create_dataset("donor_id", data=np.array(donors, dtype=object))
            obs_group.create_dataset("cell_type", data=np.array(cell_types, dtype=object))
            obs_group.create_dataset("condition_id", data=np.array(conditions, dtype=object))

            # var group
            var_group = f.create_group("var")
            var_group.attrs["_index"] = "gene_id"
            var_group.create_dataset("gene_id", data=np.array(gene_ids, dtype=object))
            var_group.create_dataset("gene_symbol", data=np.array(gene_symbols, dtype=object))

        return h5_path

    def _create_synthetic_10x_directory(
        self,
        n_cells: int = 40,
        n_genes: int = 15,
    ) -> Path:
        """Create a valid 10x Genomics Matrix Market directory."""
        dir_10x = self.dir_path / "10x_dataset"
        dir_10x.mkdir(parents=True, exist_ok=True)

        gene_ids = [f"ENSG{200000 + i:011d}" for i in range(n_genes)]
        gene_symbols = [f"FEAT_{i}" for i in range(n_genes)]
        gene_symbols[0] = "MT-ATP6"

        barcodes = [f"AAACCTGAG{j:04d}-1" for j in range(n_cells)]

        # Barcodes TSV
        with gzip.open(dir_10x / "barcodes.tsv.gz", "wt", encoding="utf-8") as f:
            for bc in barcodes:
                f.write(f"{bc}\n")

        # Features TSV
        with gzip.open(dir_10x / "features.tsv.gz", "wt", encoding="utf-8") as f:
            for gid, sym in zip(gene_ids, gene_symbols):
                f.write(f"{gid}\t{sym}\tGene Expression\n")

        # Matrix MTX (10x: genes x cells)
        rng = np.random.default_rng(123)
        mat = rng.poisson(lam=10, size=(n_genes, n_cells)).astype(np.int64)
        sp_mat = sparse.coo_matrix(mat)
        with gzip.open(dir_10x / "matrix.mtx.gz", "wb") as f:
            io.mmwrite(f, sp_mat)

        return dir_10x

    def test_h5ad_reader_extracts_genes_x_cells_matrix(self):
        """1. H5AD reader parses sparse/dense AnnData and yields SingleCellMatrix with correct orientation."""
        h5_path = self._create_synthetic_h5ad(n_cells=30, n_genes=10, is_sparse=True)

        sc_mat, obs_df, var_df = read_h5ad(h5_path)

        self.assertEqual(10, sc_mat.number_of_genes)
        self.assertEqual(30, sc_mat.number_of_cells)
        self.assertEqual(30, len(obs_df))
        self.assertEqual(10, len(var_df))
        self.assertIn("donor_id", obs_df.columns)
        self.assertIn("cell_type", obs_df.columns)
        self.assertEqual("MT-CO1", var_df.iloc[0]["gene_symbol"])

        col = sc_mat.get_cell_column(sc_mat.cell_ids[0])
        self.assertEqual(10, len(col))
        self.assertTrue(all(isinstance(int(x), int) for x in col))

    def test_10x_matrix_market_reader(self):
        """2. 10x Matrix Market directory (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz) is parsed correctly."""
        dir_10x = self._create_synthetic_10x_directory(n_cells=25, n_genes=12)

        sc_mat, obs_df, var_df = read_10x_matrix_market(dir_10x)

        self.assertEqual(12, sc_mat.number_of_genes)
        self.assertEqual(25, sc_mat.number_of_cells)
        self.assertEqual(25, len(obs_df))
        self.assertEqual(12, len(var_df))
        self.assertEqual("MT-ATP6", var_df.iloc[0]["gene_symbol"])

    def test_rds_rejection_safety(self):
        """3. R serialized files (.rds) are rejected with clear explanation and export guidance."""
        fake_rds = self.dir_path / "seurat_object.rds"
        fake_rds.write_bytes(b"RDX2\nfake_r_binary_content")

        with self.assertRaises(UnsupportedSingleCellFormatException) as ctx:
            detect_and_reject_rds(fake_rds)

        self.assertIn(".rds", str(ctx.exception))
        self.assertIn("Seurat", str(ctx.exception))
        self.assertIn("AnnData H5AD", str(ctx.exception))

    def test_single_cell_qc_evaluation_and_filtering(self):
        """4. Cell QC computes UMI counts, detected genes, mito fraction, and filters low-quality cells."""
        h5_path = self._create_synthetic_h5ad(n_cells=40, n_genes=20, is_sparse=True)
        sc_mat, obs_df, var_df = read_h5ad(h5_path)

        config = SingleCellQCConfig(min_genes=5, min_counts=20, max_mito_fraction=0.50)
        qc_res = evaluate_single_cell_qc(sc_mat, obs_df, var_df, config=config)

        self.assertEqual(40, qc_res.total_cells)
        self.assertGreater(qc_res.passed_cells, 0)
        self.assertEqual(40, len(qc_res.cell_metrics))

        first_metric = qc_res.cell_metrics[0]
        self.assertIsInstance(first_metric, CellMetadata)
        self.assertGreater(first_metric.total_counts, 0)
        self.assertGreater(first_metric.detected_genes, 0)
        self.assertGreaterEqual(first_metric.mitochondrial_fraction, 0.0)

    def test_pseudobulk_strict_raw_integer_sum_and_donor_independence(self):
        """5. Pseudobulk aggregation strictly sums raw integer counts and enforces donor independence."""
        h5_path = self._create_synthetic_h5ad(n_cells=100, n_genes=10, is_sparse=True)
        sc_mat, obs_df, var_df = read_h5ad(h5_path)

        qc_res = evaluate_single_cell_qc(sc_mat, obs_df, var_df, config=SingleCellQCConfig(min_genes=1, min_counts=1))
        # Low threshold so all populated groups qualify
        pb_config = PseudobulkConfig(min_cells_per_sample=5)

        out_tsv = self.dir_path / "pseudobulk_counts.tsv"
        pb_dataset = aggregate_pseudobulk(
            matrix=sc_mat,
            qc_result=qc_res,
            config=pb_config,
            output_tsv_path=out_tsv,
        )

        self.assertEqual(QCStatus.PASS, pb_dataset.status)
        self.assertGreater(len(pb_dataset.pseudobulk_samples), 0)
        self.assertTrue(out_tsv.exists())

        # Mathematical verification: CountMatrix values are exact integer sums of single cells
        cm = pb_dataset.count_matrix
        self.assertEqual(sc_mat.number_of_genes, len(cm.gene_ids))
        self.assertEqual(len(pb_dataset.pseudobulk_samples), len(cm.sample_ids))

        # Check a specific pseudobulk sample
        sample_0 = pb_dataset.pseudobulk_samples[0]
        # Find all cells belonging to this group
        group_cells = [
            c.cell_id for c in qc_res.cell_metrics
            if c.is_pass_qc and c.donor_id == sample_0.donor_id and c.cell_type == sample_0.cell_type and c.condition_id == sample_0.condition_id
        ]
        self.assertEqual(len(group_cells), sample_0.cell_count)

        # Direct mathematical sum check for first gene
        expected_gene0_sum = sum(sc_mat.get_value(cm.gene_ids[0], c_id) for c_id in group_cells)
        actual_gene0_sum = cm.columns[0][0]
        self.assertEqual(expected_gene0_sum, actual_gene0_sum)

        # Integer purity: all values in CountMatrix must be Python ints
        for col in cm.columns:
            for val in col:
                self.assertIsInstance(val, int)

    def test_pseudobulk_output_satisfies_frozen_bulk_intake_contract(self):
        """6. Generated pseudobulk TSV matrix is 100% compliant with Stage 6A bulk count intake validator."""
        h5_path = self._create_synthetic_h5ad(n_cells=80, n_genes=10, is_sparse=True)
        sc_mat, obs_df, var_df = read_h5ad(h5_path)

        qc_res = evaluate_single_cell_qc(sc_mat, obs_df, var_df, config=SingleCellQCConfig(min_genes=1, min_counts=1))
        out_tsv = self.dir_path / "pseudobulk_for_step6a.tsv"

        pb_dataset = aggregate_pseudobulk(
            matrix=sc_mat,
            qc_result=qc_res,
            config=PseudobulkConfig(min_cells_per_sample=2),
            output_tsv_path=out_tsv,
        )

        # Validate with frozen Stage 6A matrix reader
        read_cm, findings = read_matrix(out_tsv, {})
        self.assertIsNotNone(read_cm)
        self.assertEqual(len(pb_dataset.count_matrix.sample_ids), len(read_cm.sample_ids))

        # Validate with Stage 6A validate_bulk_counts
        qc_result_6a = validate_bulk_counts(
            path=out_tsv,
            metadata=pb_dataset.metadata,
            asset_id=pb_dataset.metadata["input_data_assets"][0]["asset_id"],
            assay_id=pb_dataset.metadata["sequencing_assays"][0]["assay_id"],
        )
        self.assertEqual(QCStatus.PASS, qc_result_6a.status)

    def test_dispatcher_planning_and_execution_end_to_end(self):
        """7. Dispatcher activates single_cell_pipeline_runner and hands off pseudobulk to bulk pipeline."""
        h5_path = self._create_synthetic_h5ad(n_cells=60, n_genes=10, is_sparse=True)

        meta = routing_metadata(h5_path, "H5AD", modality="scRNA_seq")
        meta["enable_single_cell_execution"] = True
        meta["output_tsv_path"] = str(self.dir_path / "dispatcher_pseudobulk.tsv")

        plan = plan_route_execution(meta, enable_u4_single_cell_execution=True)

        self.assertEqual(Route.SINGLE_CELL_OBJECT_PIPELINE, plan.route)
        self.assertEqual("READY_FOR_EXECUTION", plan.status)
        self.assertEqual("single_cell_pipeline_runner", plan.executor_id)
        self.assertEqual("pseudobulk_raw_counts", plan.target_format)
        self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, plan.downstream_route)

        res = dispatch_route_execution(
            plan,
            qc_config=SingleCellQCConfig(min_genes=1, min_counts=1),
            pseudobulk_config=PseudobulkConfig(min_cells_per_sample=2),
            output_tsv_path=meta["output_tsv_path"],
        )

        self.assertTrue(res.is_success)
        self.assertEqual("COMPLETED", res.execution_state)
        self.assertEqual("pseudobulk_dataset", res.emitted_format)
        self.assertIn("dataset", res.artifacts)
        self.assertIn("count_matrix", res.artifacts)

        pb_ds = res.artifacts["dataset"]
        self.assertIsInstance(pb_ds, PseudobulkDataset)
        self.assertGreater(len(pb_ds.pseudobulk_samples), 0)


if __name__ == "__main__":
    unittest.main()
