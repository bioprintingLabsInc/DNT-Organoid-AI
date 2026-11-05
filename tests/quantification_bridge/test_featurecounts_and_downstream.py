"""Unit and integration tests for featureCounts runner and downstream bridge.

Verifies:
1. featureCounts command construction (-p --countReadPairs -B -C -s, omitting -M and -O).
2. Summary metric parsing from .summary output.
3. Count matrix adapter conversion and numerical conservation.
4. Canonical metadata generation for PRJNA802072.
5. Biological replicate guard in downstream pipeline (blocks DESeq2 when n < 2).
"""

import json
from pathlib import Path
import tempfile
import pytest

from src.quantification_bridge.downstream_bridge import (
    DownstreamBridgeResult,
    build_canonical_metadata_for_runs,
    execute_downstream_pipeline,
)
from src.quantification_bridge.featurecounts_runner import (
    parse_featurecounts_summary,
    run_featurecounts,
)
from src.quantification_bridge.adapter import convert_featurecounts_to_matrix
from src.quantification_bridge.manifest import (
    APPROVED_2_SAMPLE_TECHNICAL_PILOT,
    APPROVED_PILOT_RUNS,
)


def test_parse_featurecounts_summary(tmp_path: Path):
    """Verify that featureCounts .summary files are accurately parsed into metrics."""
    summary_file = tmp_path / "counts.txt.summary"
    content = (
        "Status\tsample1_Aligned.sortedByCoord.out.bam\tsample2_Aligned.sortedByCoord.out.bam\n"
        "Assigned\t18500200\t19200150\n"
        "Unassigned_Unmapped\t500000\t450000\n"
        "Unassigned_Read_Type\t0\t0\n"
        "Unassigned_Singleton\t0\t0\n"
        "Unassigned_MappingQuality\t0\t0\n"
        "Unassigned_Chimera\t12000\t15000\n"
        "Unassigned_FragmentLength\t0\t0\n"
        "Unassigned_Duplicate\t0\t0\n"
        "Unassigned_MultiMapping\t1450000\t1380000\n"
        "Unassigned_Secondary\t0\t0\n"
        "Unassigned_NonSplit\t0\t0\n"
        "Unassigned_NoFeatures\t850000\t910000\n"
        "Unassigned_Overlapping_Length\t0\t0\n"
        "Unassigned_Ambiguity\t210000\t225000\n"
    )
    summary_file.write_text(content, encoding="utf-8")

    parsed = parse_featurecounts_summary(summary_file)
    assert len(parsed) == 2
    s1 = "sample1_Aligned.sortedByCoord.out.bam"
    s2 = "sample2_Aligned.sortedByCoord.out.bam"
    assert s1 in parsed and s2 in parsed
    assert parsed[s1]["Assigned"] == 18500200
    assert parsed[s1]["Unassigned_MultiMapping"] == 1450000
    assert parsed[s1]["Unassigned_NoFeatures"] == 850000
    assert parsed[s1]["Unassigned_Ambiguity"] == 210000

    assert parsed[s2]["Assigned"] == 19200150
    assert parsed[s2]["Unassigned_MultiMapping"] == 1380000


def test_featurecounts_invalid_strandedness():
    """Verify that invalid strandedness flags raise a ValueError."""
    with pytest.raises(ValueError, match="Invalid strandedness flag"):
        run_featurecounts(
            bam_paths=[],
            gtf_path=Path("dummy.gtf"),
            output_counts_path=Path("out.txt"),
            strandedness=3,
        )


def test_convert_featurecounts_to_matrix_synthetic(tmp_path: Path):
    """Verify adapter conversion preserves 100% numerical values without modification."""
    fc_out = tmp_path / "raw_featurecounts.txt"
    matrix_out = tmp_path / "matrix.tsv"
    universe_out = tmp_path / "universe.json"

    # Synthetic featureCounts output
    lines = [
        "# Program:featureCounts v2.0.6; Command:featureCounts ...",
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/data/align/SRR17834920_Aligned.sortedByCoord.out.bam\t/data/align/SRR17834905_Aligned.sortedByCoord.out.bam",
        "ENSG00000000003\tchrX\t100627108\t100639991\t-\t4535\t142\t89",
        "ENSG00000000005\tchrX\t100584802\t100599885\t+\t1476\t0\t0",
        "ENSG00000000419\tchr20\t50934867\t50958555\t-\t1207\t3560\t3120",
        "ENSG00000000457\tchr1\t169818772\t169863408\t-\t6361\t890\t940",
    ]
    fc_out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    approved_runs = ("SRR17834920", "SRR17834905")
    res = convert_featurecounts_to_matrix(
        raw_featurecounts_path=fc_out,
        output_matrix_path=matrix_out,
        output_universe_path=universe_out,
        approved_runs=approved_runs,
    )

    assert res["is_numerically_conserved"] is True
    assert res["total_genes_emitted"] == 4
    assert res["sample_columns_emitted"] == ["SRR17834920", "SRR17834905"]

    matrix_lines = matrix_out.read_text(encoding="utf-8").strip().splitlines()
    assert matrix_lines[0] == "gene_id\tSRR17834920\tSRR17834905"
    assert matrix_lines[1] == "ENSG00000000003\t142\t89"
    assert matrix_lines[2] == "ENSG00000000005\t0\t0"
    assert matrix_lines[3] == "ENSG00000000419\t3560\t3120"
    assert matrix_lines[4] == "ENSG00000000457\t890\t940"

    universe_lines = universe_out.read_text(encoding="utf-8").strip().splitlines()
    assert universe_lines[0] == "gene_id\tchr\tstart\tend\tstrand\tlength\ttotal_counts"
    assert len(universe_lines) == 5
    assert "ENSG00000000003" in universe_lines[1]


def test_build_canonical_metadata_for_runs():
    """Verify canonical metadata contains all required biological and assay fields."""
    runs = ("SRR17834920", "SRR17834905")
    meta = build_canonical_metadata_for_runs(runs)

    assert meta["studies"][0]["study_id"] == "PRJNA802072"
    assert len(meta["samples"]) == 2
    assert len(meta["sequencing_assays"]) == 1

    # Check Control sample (SRR17834920)
    ctrl = next(s for s in meta["samples"] if "SRR17834920" in s["sample_id"])
    assert ctrl["treatment_group"] == "control"
    assert ctrl["compound_concentration_um"]["original_value"] == 0.0
    assert ctrl["developmental_day"]["original_value"] == 28
    assert ctrl["tissue_or_cell_type"]["original_value"] == "Human Brain Cortical Organoid"

    # Check Treatment sample (SRR17834905)
    trt = next(s for s in meta["samples"] if "SRR17834905" in s["sample_id"])
    assert trt["treatment_group"] == "treatment"
    assert trt["compound_concentration_um"]["original_value"] == 250.0
    assert trt["developmental_day"]["original_value"] == 28


def test_downstream_biological_replication_guard(tmp_path: Path):
    """Verify that downstream pipeline blocks DESeq2 when biological replicates n < 2."""
    matrix_file = tmp_path / "counts_pilot.tsv"
    # Create minimal count matrix for 2 samples
    lines = [
        "gene_id\tSRR17834920\tSRR17834905",
        "ENSG00000000003\t100\t120",
        "ENSG00000000419\t500\t480",
        "ENSG00000000457\t250\t300",
    ]
    matrix_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=("SRR17834920", "SRR17834905"),
        output_dir=tmp_path / "downstream",
    )

    # Harmonization and QC should have evaluated
    assert result.status in ("DOWNSTREAM_PARTIAL_INSUFFICIENT_REPLICATION", "DOWNSTREAM_EVALUATION_COMPLETED")
    assert result.qc_passed is True

    # DESeq2 MUST be blocked by biological replicate guard
    assert result.deseq2_status == "BLOCKED_INSUFFICIENT_BIOLOGICAL_REPLICATION"
    assert result.biological_replicate_warning is not None
    assert "n >= 2 biological replicates" in result.biological_replicate_warning
    assert result.pathway_status == "NOT_EVALUATED"


def test_downstream_bridge_metadata_validation_and_integrity_gates(tmp_path: Path):
    """Verify that MetadataValidator and MetadataIntegrityChecker run during downstream bridge execution."""
    matrix_file = tmp_path / "counts_pilot.tsv"
    lines = [
        "gene_id\tSRR17834920\tSRR17834905",
        "ENSG00000000003\t100\t120",
        "ENSG00000000419\t500\t480",
        "ENSG00000000457\t250\t300",
    ]
    matrix_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=("SRR17834920", "SRR17834905"),
        output_dir=tmp_path / "downstream",
    )

    # Metadata validator must have executed
    assert result.metadata_validation_status.lower() in ("valid_with_warnings", "valid")
    assert result.metadata_validation_metrics is not None
    assert result.metadata_validation_metrics["validator_status"].lower() in ("valid_with_warnings", "valid")
    assert result.metadata_validation_metrics["integrity_is_valid"] is False
    assert result.metadata_validation_metrics["integrity_status"].lower() == "invalid"
    assert result.metadata_validation_metrics["integrity_findings_count"] > 0
    assert "replicate_summary" in result.metadata_validation_metrics


def _generate_synthetic_eligible_count_matrix(matrix_path: Path, runs: tuple[str, ...], n_genes: int = 100) -> None:
    """Helper to generate an independently replicated count matrix using real Ensembl/Entrez genes."""
    import random
    from src.features.pathway.mapping import EnsemblToEntrezAdapter
    from src.quantification_bridge.manifest import RUN_TO_SAMPLE_MAP

    adapter = EnsemblToEntrezAdapter()
    gene_recs = [
        r for r in adapter.gene_reference.records
        if r.entrez_gene_id and "|" not in r.entrez_gene_id and r.gene_biotype == "protein_coding"
    ][:n_genes]

    rng = random.Random(42)
    header = "gene_id\t" + "\t".join(runs)
    lines = [header]
    ctrl_runs = [r for r in runs if RUN_TO_SAMPLE_MAP.get(r) and RUN_TO_SAMPLE_MAP[r].condition_group == "control"]
    trt_runs = [r for r in runs if RUN_TO_SAMPLE_MAP.get(r) and RUN_TO_SAMPLE_MAP[r].condition_group == "treatment"]

    for i, r in enumerate(gene_recs):
        base = 100 + (i * 15)
        fold = 2.5 if i < 25 else (0.4 if i < 50 else 1.0)
        row = [r.ensembl_gene_id]
        for run_id in runs:
            is_trt = run_id in trt_runs
            m = fold if is_trt else 1.0
            count = max(10, int(base * m * rng.uniform(0.9, 1.1)))
            row.append(str(count))
        lines.append("\t".join(row))

    matrix_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_downstream_bridge_verified_biological_replicates_eligible_outcome(tmp_path: Path):
    """Verify that an eligible cohort with verified biological replicates completes DESeq2 and Pathway analysis."""
    matrix_file = tmp_path / "counts_verified.tsv"
    runs = ("SRR17834920", "SRR17834908", "SRR17834905", "SRR17834904")
    _generate_synthetic_eligible_count_matrix(matrix_file, runs, n_genes=100)

    verified_reps = {
        "SRR17834920": "c_rep1",
        "SRR17834908": "c_rep2",
        "SRR17834905": "t_rep1",
        "SRR17834904": "t_rep2",
    }

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        output_dir=tmp_path / "downstream_verified",
        verified_biological_replicates=verified_reps,
    )

    # 1. Exploratory QC and normalization pass
    assert result.qc_passed is True
    assert result.normalization_metrics["status"] in ("PASS", "Status.PASS")

    # 2. Integrity checker verifies biological replication
    assert result.metadata_validation_metrics["integrity_is_valid"] is True
    assert result.metadata_validation_metrics["integrity_status"].lower() == "valid"

    # 3. Explicitly verify DESEQ2_COMPLETED (strengthened check)
    assert result.deseq2_status == "DESEQ2_COMPLETED"
    assert result.deseq2_metrics["status"] == "DESEQ2_SUCCESS"
    assert result.deseq2_metrics["tested_genes"] == 100
    assert result.deseq2_metrics["significant_degs"] > 0
    assert result.deseq2_metrics["design_formula"] == "~ condition"

    # 4. Explicitly verify Stage 8C Pathway Analysis output contract
    assert result.pathway_status == "PATHWAY_ANALYSIS_SUCCESS"
    assert result.pathway_metrics["status"] == "GSEA_SUCCESS"
    assert result.pathway_metrics["qc_status"] in ("PASS", "WARNING")
    assert result.pathway_metrics["ranked_genes_count"] == 100
    assert result.pathway_metrics["total_eligible_genes"] == 100
    assert result.pathway_metrics["go_results_count"] > 0
    assert "significant_pathways" in result.pathway_metrics

    # 5. Full workflow status distinguishes completed analysis
    assert result.status == "DOWNSTREAM_COMPLETED"


def test_downstream_bridge_invalid_canonical_metadata_blocks_inferential_analysis(tmp_path: Path):
    """Verify that invalid canonical metadata blocks inferential analysis while permitting exploratory QC/normalization."""
    matrix_file = tmp_path / "counts_invalid_meta.tsv"
    runs = ("SRR17834920", "SRR17834908", "SRR17834905", "SRR17834904")
    _generate_synthetic_eligible_count_matrix(matrix_file, runs, n_genes=50)

    # Construct invalid canonical metadata (empty studies list violates schema)
    invalid_meta = build_canonical_metadata_for_runs(runs, strict_schema=True)
    invalid_meta["studies"] = []

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        output_dir=tmp_path / "downstream_invalid_meta",
        canonical_metadata=invalid_meta,
    )

    # Exploratory processing permitted on count matrix
    assert result.qc_passed is True

    # Metadata validator detects invalid metadata
    assert result.metadata_validation_status.lower() == "invalid"
    assert result.metadata_validation_metrics["validator_is_valid"] is False

    # Inferential analysis is blocked and reasons recorded
    assert result.status == "DOWNSTREAM_BLOCKED_INVALID_METADATA"
    assert result.deseq2_status == "BLOCKED_INVALID_CANONICAL_METADATA"
    assert result.pathway_status == "BLOCKED_INVALID_CANONICAL_METADATA"
    assert result.biological_replicate_warning is not None
    assert "Canonical metadata validation failed" in result.biological_replicate_warning
    assert "findings" in result.deseq2_metrics
    assert len(result.deseq2_metrics["findings"]) > 0


def test_downstream_bridge_pathway_output_contract_on_qc_failure(tmp_path: Path):
    """Verify pathway analysis output contract does not accept placeholder values or success when QC fails."""
    matrix_file = tmp_path / "counts_nonentrez.tsv"
    runs = ("SRR17834920", "SRR17834908", "SRR17834905", "SRR17834904")

    # Generate 100 genuine Ensembl 112 genes that lack Entrez IDs (so 0 Entrez mappings exist)
    from src.features.pathway.mapping import EnsemblToEntrezAdapter
    adapter = EnsemblToEntrezAdapter()
    gene_recs = [r for r in adapter.gene_reference.records if not r.entrez_gene_id][:100]

    import random
    rng = random.Random(99)
    lines = ["gene_id\t" + "\t".join(runs)]
    for i, r in enumerate(gene_recs):
        base = 100 + i * 15
        fold = 2.0 if i < 20 else (0.5 if i < 40 else 1.0)
        c1 = max(10, int(base * rng.uniform(0.9, 1.1)))
        c2 = max(10, int(base * rng.uniform(0.9, 1.1)))
        t1 = max(10, int(base * fold * rng.uniform(0.9, 1.1)))
        t2 = max(10, int(base * fold * rng.uniform(0.9, 1.1)))
        lines.append(f"{r.ensembl_gene_id}\t{c1}\t{c2}\t{t1}\t{t2}")
    matrix_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    verified_reps = {
        "SRR17834920": "c_rep1",
        "SRR17834908": "c_rep2",
        "SRR17834905": "t_rep1",
        "SRR17834904": "t_rep2",
    }

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        output_dir=tmp_path / "downstream_qc_fail",
        verified_biological_replicates=verified_reps,
    )

    # QC and normalization pass
    assert result.qc_passed is True

    # DESeq2 completes successfully
    assert result.deseq2_status == "DESEQ2_COMPLETED"

    # Pathway analysis evaluated: 0 genes mapped to Entrez, so QC must FAIL
    # It must NOT return PATHWAY_ANALYSIS_SUCCESS or placeholders
    assert result.pathway_status == "PATHWAY_ANALYSIS_QC_FAIL"
    assert result.pathway_metrics["qc_status"] == "FAIL"
    assert result.pathway_metrics["status"] == "GSEA_FAIL"
    assert result.pathway_metrics["ranked_genes_count"] == 0
    assert result.pathway_metrics["go_results_count"] == 0
    assert result.pathway_metrics["kegg_results_count"] == 0
    assert len(result.pathway_metrics["qc_notes"]) > 0

    # Workflow status distinguishes partially completed analysis (DESeq2 ok, pathway failed)
    assert result.status == "DOWNSTREAM_PARTIAL_PATHWAY_FAILED"
