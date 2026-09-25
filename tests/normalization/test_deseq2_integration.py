"""Live scientific validation and integration tests using official Bioconductor DESeq2.

Exercised against the live installed R 4.4.3, Bioconductor 3.20, DESeq2 1.46.0 runtime.
Zero Python approximations used.
"""

from typing import Any
import pytest

from src.gene_harmonization.models import (
    GeneMappingStatus,
    HarmonizationSummary,
    HarmonizedDataset,
    HarmonizedGene,
    Status as HarmonizationStatus,
)
from src.normalization.config import NormalizationConfig
from src.normalization.models import (
    NormalizationMethod,
    NormalizationReviewDisposition,
    ResponseEligibility,
    Severity,
    Status,
)
from src.normalization.normalizer import normalize_bulk_dataset
from src.normalization.r_bridge import probe_r_environment
from src.qc.matrix import CountMatrix
from src.qc.result import Status as QCStatus


def test_r_environment_locked_versions() -> None:
    """Validate that the active R environment reports exactly the version-locked baseline."""
    env_info = probe_r_environment("Rscript")
    assert env_info.is_available, f"DESeq2 environment is unavailable: {env_info.error_message}"
    assert env_info.r_version == "4.4.3", f"Expected R 4.4.3, got {env_info.r_version}"
    assert env_info.bioc_version == "3.20", f"Expected Bioconductor 3.20, got {env_info.bioc_version}"
    assert env_info.deseq2_version == "1.46.0", f"Expected DESeq2 1.46.0, got {env_info.deseq2_version}"


def create_benchmark_fixture() -> tuple[HarmonizedDataset, dict[str, Any]]:
    """Create fixed benchmark integration fixture: 5 genes, 3 samples.

    Sample 2 has exactly 2.0x the sequencing depth of Sample 1 and Sample 3.
    """
    sample_ids = ("sample_1", "sample_2", "sample_3")
    col1 = (100, 200, 300, 400, 500)
    col2 = (200, 400, 600, 800, 1000)
    col3 = (100, 200, 300, 400, 500)
    columns = (col1, col2, col3)
    gene_ids = ("ENSG00000000001", "ENSG00000000002", "ENSG00000000003", "ENSG00000000004", "ENSG00000000005")

    matrix = CountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=columns)

    h_genes = tuple(
        HarmonizedGene(
            source_index=i,
            original_gene_id=gid,
            identifier_type="ensembl_gene_id",
            canonical_gene_id=gid,
            approved_symbol=f"GENE_{i+1}",
            mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
            mapping_reason=None,
            candidate_canonical_ids=(),
            reference_id="ensembl_human_genes_v112",
            reference_version="112",
            source_asset_id="asset_benchmark_001",
        )
        for i, gid in enumerate(gene_ids)
    )

    summary = HarmonizationSummary(
        total_genes=5,
        uniquely_mapped_genes=5,
        unmapped_genes=0,
        ambiguous_genes=0,
        invalid_genes=0,
        canonical_collision_genes=0,
        colliding_source_gene_count=0,
        percentage_uniquely_mapped=100.0,
        percentage_unresolved=0.0,
    )

    harm = HarmonizedDataset(
        status=HarmonizationStatus.PASS,
        findings=(),
        summary=summary,
        genes=h_genes,
        collisions=(),
        matrix=matrix,
        source_asset_id="asset_benchmark_001",
        assay_id="assay_001",
        sample_ids=sample_ids,
        reference_identity="Ensembl Release 112",
        reference_version="112",
        reference_checksum="3054aa734c15bb0a9cc64b67a8ead2a023e7dfd2cc0d0d6b7706372f4b13940e",
        config_version="1.0.0",
        config_checksum="abc",
        harmonization_version="bulk_gene_harmonization_v1",
        upstream_qc_status=QCStatus.PASS,
    )

    meta = {
        "studies": [{"study_id": "study_bench"}],
        "experiments": [{"experiment_id": "exp_bench", "study_id": "study_bench"}],
        "conditions": [
            {"condition_id": "cond_treat", "experiment_id": "exp_bench", "treatment_control_status": "treatment"},
            {"condition_id": "cond_ctrl", "experiment_id": "exp_bench", "treatment_control_status": "control"},
        ],
        "samples": [
            {"sample_id": "sample_1", "experiment_id": "exp_bench", "condition_id": "cond_treat"},
            {"sample_id": "sample_2", "experiment_id": "exp_bench", "condition_id": "cond_treat"},
            {"sample_id": "sample_3", "experiment_id": "exp_bench", "condition_id": "cond_ctrl"},
        ],
        "sequencing_assays": [
            {
                "assay_id": "assay_001",
                "sample_ids": ["sample_1", "sample_2", "sample_3"],
                "rna_seq_modality": "bulk_rna_seq",
                "library_strategy": "polyA",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "strandedness": "reverse",
                "read_layout": "paired_end",
            }
        ],
        "treatment_control_relationships": [
            {
                "relationship_id": "rel_bench",
                "treatment_condition_ids": ["cond_treat"],
                "matched_control_condition_ids": ["cond_ctrl"],
            }
        ],
    }

    return harm, meta


def test_live_official_deseq2_execution_and_numerical_parity() -> None:
    """Live validation of official DESeq2 size factor estimation and numerical parity."""
    harm, meta = create_benchmark_fixture()
    raw_col1 = harm.matrix.columns[0]
    raw_col2 = harm.matrix.columns[1]
    raw_col3 = harm.matrix.columns[2]

    res = normalize_bulk_dataset(harm, meta)

    assert res.status == Status.PASS
    assert len(res.cohort_results) == 1
    cohort_res = res.cohort_results[0]

    # 1. Verify input raw integer counts are unchanged
    assert cohort_res.raw_count_matrix.columns[0] == raw_col1
    assert cohort_res.raw_count_matrix.columns[1] == raw_col2
    assert cohort_res.raw_count_matrix.columns[2] == raw_col3

    # 2. Verify official DESeq2 size factors
    sf_dict = {sf.sample_id: sf.size_factor for sf in cohort_res.size_factors}
    assert len(sf_dict) == 3
    assert all(sf > 0 for sf in sf_dict.values())

    # Analytical calculation:
    # geo_means = (100 * 200 * 100)^(1/3) = 125.992104989
    # sample_1: 100 / 125.992104989 = 0.7937005
    # sample_2: 200 / 125.992104989 = 1.5874011
    # sample_3: 100 / 125.992104989 = 0.7937005
    assert pytest.approx(sf_dict["sample_1"], rel=1e-5) == 0.7937005
    assert pytest.approx(sf_dict["sample_2"], rel=1e-5) == 1.5874011
    assert pytest.approx(sf_dict["sample_3"], rel=1e-5) == 0.7937005
    assert pytest.approx(sf_dict["sample_2"] / sf_dict["sample_1"], rel=1e-5) == 2.0

    # 3. Verify diagnostic normalized matrix: K_ij / s_j
    norm_mat = cohort_res.diagnostic_normalized_matrix
    assert norm_mat is not None
    for j, sid in enumerate(cohort_res.cohort.sample_ids):
        sf = sf_dict[sid]
        for i in range(len(norm_mat.gene_ids)):
            expected_norm = float(cohort_res.raw_count_matrix.columns[j][i]) / sf
            actual_norm = norm_mat.columns[j][i]
            assert pytest.approx(actual_norm, rel=1e-5) == expected_norm

    # 4. Verify treatment/control membership & response eligibility
    sf_objs = {sf.sample_id: sf for sf in cohort_res.size_factors}
    assert sf_objs["sample_1"].treatment_control_status == "treatment"
    assert sf_objs["sample_1"].response_eligibility == ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL
    assert sf_objs["sample_1"].matched_control_sample_ids == ("sample_3",)

    assert sf_objs["sample_2"].treatment_control_status == "treatment"
    assert sf_objs["sample_2"].response_eligibility == ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL
    assert sf_objs["sample_2"].matched_control_sample_ids == ("sample_3",)

    assert sf_objs["sample_3"].treatment_control_status == "control"
    assert sf_objs["sample_3"].response_eligibility == ResponseEligibility.CONTROL_BASELINE

    # 5. Verify provenance contains exact versions
    env_info = cohort_res.r_environment_info
    assert env_info["r_version"] == "4.4.3"
    assert env_info["bioc_version"] == "3.20"
    assert env_info["deseq2_version"] == "1.46.0"
    assert env_info["is_available"] is True

    # 6. Verify numerical repeatability
    res_repeat = normalize_bulk_dataset(harm, meta)
    sf_repeat = {sf.sample_id: sf.size_factor for sf in res_repeat.cohort_results[0].size_factors}
    for sid in sf_dict:
        assert pytest.approx(sf_dict[sid], rel=1e-10) == sf_repeat[sid]


def test_live_standard_ratio_zero_geomean_failure_and_poscounts_gating() -> None:
    """Exercise live DESeq2 where standard ratio cannot compute size factors, and verify poscounts authorization."""
    # Matrix where every gene has at least one zero across samples
    # Gene 1: (0, 100, 200)
    # Gene 2: (100, 0, 200)
    # Gene 3: (100, 200, 0)
    sample_ids = ("sample_1", "sample_2", "sample_3")
    col1 = (0, 100, 100)
    col2 = (100, 0, 200)
    col3 = (200, 200, 0)
    columns = (col1, col2, col3)
    gene_ids = ("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")

    matrix = CountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=columns)

    h_genes = tuple(
        HarmonizedGene(
            source_index=i,
            original_gene_id=gid,
            identifier_type="ensembl_gene_id",
            canonical_gene_id=gid,
            approved_symbol=f"GENE_{i+1}",
            mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
            mapping_reason=None,
            candidate_canonical_ids=(),
            reference_id="ensembl_human_genes_v112",
            reference_version="112",
            source_asset_id="asset_sparse",
        )
        for i, gid in enumerate(gene_ids)
    )

    summary = HarmonizationSummary(
        total_genes=3,
        uniquely_mapped_genes=3,
        unmapped_genes=0,
        ambiguous_genes=0,
        invalid_genes=0,
        canonical_collision_genes=0,
        colliding_source_gene_count=0,
        percentage_uniquely_mapped=100.0,
        percentage_unresolved=0.0,
    )

    harm = HarmonizedDataset(
        status=HarmonizationStatus.PASS,
        findings=(),
        summary=summary,
        genes=h_genes,
        collisions=(),
        matrix=matrix,
        source_asset_id="asset_sparse",
        assay_id="assay_sparse",
        sample_ids=sample_ids,
        reference_identity="Ensembl Release 112",
        reference_version="112",
        reference_checksum="abc",
        config_version="1.0.0",
        config_checksum="def",
        harmonization_version="bulk_gene_harmonization_v1",
        upstream_qc_status=QCStatus.PASS,
    )

    meta = {
        "studies": [{"study_id": "study_sparse"}],
        "experiments": [{"experiment_id": "exp_sparse", "study_id": "study_sparse"}],
        "conditions": [
            {"condition_id": "cond_treat", "experiment_id": "exp_sparse", "treatment_control_status": "treatment"},
            {"condition_id": "cond_ctrl", "experiment_id": "exp_sparse", "treatment_control_status": "control"},
        ],
        "samples": [
            {"sample_id": "sample_1", "experiment_id": "exp_sparse", "condition_id": "cond_treat"},
            {"sample_id": "sample_2", "experiment_id": "exp_sparse", "condition_id": "cond_treat"},
            {"sample_id": "sample_3", "experiment_id": "exp_sparse", "condition_id": "cond_ctrl"},
        ],
        "sequencing_assays": [
            {
                "assay_id": "assay_sparse",
                "sample_ids": ["sample_1", "sample_2", "sample_3"],
                "rna_seq_modality": "bulk_rna_seq",
                "library_strategy": "polyA",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "strandedness": "reverse",
                "read_layout": "paired_end",
            }
        ],
        "treatment_control_relationships": [
            {
                "relationship_id": "rel_sparse",
                "treatment_condition_ids": ["cond_treat"],
                "matched_control_condition_ids": ["cond_ctrl"],
            }
        ],
    }

    # 1. Run with default ratio method: MUST FAIL WITH REVIEW FINDING, NO SILENT FALLBACK
    res_ratio = normalize_bulk_dataset(harm, meta)

    assert res_ratio.status == Status.NEEDS_REVIEW
    assert any(f.rule_id == "standard_size_factors_failed" for f in res_ratio.findings)
    assert any("every gene contains at least one zero" in f.message.lower() for f in res_ratio.findings)

    cohort_res = res_ratio.cohort_results[0]
    assert cohort_res.status == Status.NEEDS_REVIEW
    assert cohort_res.diagnostic_normalized_matrix is None
    # Raw counts remain available and intact
    assert cohort_res.raw_count_matrix.columns == columns

    # 2. Run with explicitly authorized poscounts
    pos_disp = NormalizationReviewDisposition(
        decision="AUTHORIZE_POSCOUNTS",
        rationale="Sparse benchmark matrix where every gene contains zero across samples; poscounts authorized.",
        reviewer="scientific_curator",
        review_date="2026-09-24",
        affected_cohort_id=cohort_res.cohort.cohort_id,
    )
    cfg_pos = NormalizationConfig(size_factor_method="poscounts")

    res_pos = normalize_bulk_dataset(harm, meta, config=cfg_pos, review_disposition=pos_disp)

    assert res_pos.status == Status.PASS
    pos_cohort_res = res_pos.cohort_results[0]
    assert pos_cohort_res.status == Status.PASS
    assert pos_cohort_res.method == "poscounts"
    assert pos_cohort_res.diagnostic_normalized_matrix is not None
    assert len(pos_cohort_res.size_factors) == 3
    assert all(sf.size_factor > 0 for sf in pos_cohort_res.size_factors)


def test_live_normalization_cohort_preservation_and_isolation() -> None:
    """Confirm live DESeq2 execution respects cohort boundaries and preserves matched treatment/controls."""
    # Treatment is on lane_1, control is on lane_2
    harm, _ = create_benchmark_fixture()

    meta_multi_assay = {
        "studies": [{"study_id": "study_bench"}],
        "experiments": [{"experiment_id": "exp_bench", "study_id": "study_bench"}],
        "conditions": [
            {"condition_id": "cond_treat", "experiment_id": "exp_bench", "treatment_control_status": "treatment"},
            {"condition_id": "cond_ctrl", "experiment_id": "exp_bench", "treatment_control_status": "control"},
        ],
        "samples": [
            {"sample_id": "sample_1", "experiment_id": "exp_bench", "condition_id": "cond_treat"},
            {"sample_id": "sample_2", "experiment_id": "exp_bench", "condition_id": "cond_treat"},
            {"sample_id": "sample_3", "experiment_id": "exp_bench", "condition_id": "cond_ctrl"},
        ],
        "sequencing_assays": [
            {
                "assay_id": "assay_lane_1",
                "sample_ids": ["sample_1", "sample_2"],
                "rna_seq_modality": "bulk_rna_seq",
                "library_strategy": "polyA",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "strandedness": "reverse",
                "read_layout": "paired_end",
            },
            {
                "assay_id": "assay_lane_2",
                "sample_ids": ["sample_3"],
                "rna_seq_modality": "bulk_rna_seq",
                "library_strategy": "polyA",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "strandedness": "reverse",
                "read_layout": "paired_end",
            },
        ],
        "treatment_control_relationships": [
            {
                "relationship_id": "rel_bench",
                "treatment_condition_ids": ["cond_treat"],
                "matched_control_condition_ids": ["cond_ctrl"],
            }
        ],
    }

    res = normalize_bulk_dataset(harm, meta_multi_assay)
    assert res.status == Status.PASS
    # Treatment and matched controls spanning assays lane_1 and lane_2 must be co-normalized in the SAME cohort
    assert len(res.cohort_results) == 1
    cohort = res.cohort_results[0].cohort
    assert set(cohort.sample_ids) == {"sample_1", "sample_2", "sample_3"}
    assert set(cohort.assay_ids) == {"assay_lane_1", "assay_lane_2"}


def test_missing_r_environment_reporting_cleanly() -> None:
    """Verify that when an invalid R executable path is supplied, the pipeline reports r_environment_missing without fake Python fallbacks."""
    harm, meta = create_benchmark_fixture()
    cfg_invalid = NormalizationConfig(r_binary_path="nonexistent_r_path_xyz_123")

    res = normalize_bulk_dataset(harm, meta, config=cfg_invalid)

    assert res.status == Status.FAIL
    assert any(f.rule_id == "r_environment_missing" for f in res.findings)
    assert any(f.severity == Severity.ERROR for f in res.findings)
    assert len(res.cohort_results) == 1
    assert res.cohort_results[0].status == Status.FAIL
    assert res.cohort_results[0].diagnostic_normalized_matrix is None
    assert harm.matrix is not None
