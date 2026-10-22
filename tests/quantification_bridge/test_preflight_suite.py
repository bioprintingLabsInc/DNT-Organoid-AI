"""Mandatory synthetic preflight validation suite for PRJNA802072 quantification bridge.

Covers tests A through G required before sequencing data execution:
A. Paired-fragment counting test (-p --countReadPairs -B -C)
B. Multimapper exclusion test (no -M, NH:i:2)
C. Ambiguous overlap exclusion test (no -O -> Unassigned_Ambiguity)
D. Adapter numerical-conservation test (0 discrepancies)
E. Sample-header mapping test
F. BED12 reproducibility test (deterministic SHA-256)
G. Frozen-pipeline immutability test (git diff check)
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path

import pytest

from src.quantification_bridge.adapter import (
    convert_featurecounts_to_matrix,
    map_bam_header_to_run_id,
    sha256_file,
)
from src.quantification_bridge.bed_generator import generate_bed12_from_gtf
from src.quantification_bridge.env_exec import run_wsl, to_wsl_path
from src.quantification_bridge.manifest import APPROVED_RUN_ACCESSIONS


def test_a_b_c_featurecounts_synthetic_bam():
    """Execute Tests A, B, and C using a minimal synthetic BAM and GTF."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        sam_path = tmp_p / "synthetic.sam"
        bam_path = tmp_p / "synthetic.bam"
        sorted_bam_path = tmp_p / "synthetic_sorted.bam"
        gtf_path = tmp_p / "synthetic.gtf"
        fc_out_path = tmp_p / "fc_counts.txt"

        # 1. Create Synthetic GTF
        # GeneA: chr1:100-500
        # GeneB: chr1:1000-1500 (overlaps GeneC)
        # GeneC: chr1:1200-1800 (overlaps GeneB)
        # GeneD: chr2:100-500
        gtf_content = (
            'chr1\tsynth\texon\t100\t500\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";\n'
            'chr1\tsynth\texon\t1000\t1500\t.\t+\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002";\n'
            'chr1\tsynth\texon\t1200\t1800\t.\t+\t.\tgene_id "ENSG00000000003"; transcript_id "ENST00000000003";\n'
            'chr2\tsynth\texon\t100\t500\t.\t+\t.\tgene_id "ENSG00000000004"; transcript_id "ENST00000000004";\n'
        )
        gtf_path.write_text(gtf_content, encoding="utf-8")

        # 2. Create Synthetic SAM
        # Header
        # Reads:
        # pair1: mapped to GeneA (pos 150-350) -> valid concordant pair (+1)
        # pair2: mapped to GeneA (pos 200-400) -> valid concordant pair (+1)
        # orphan1: R1 at pos 250 (GeneA), R2 unmapped -> rejected by -B
        # chimeric1: R1 at chr1:220 (GeneA), R2 at chr2:200 (GeneD) -> rejected by -C
        # multimap1: R1 at pos 260, R2 at pos 400, NH:i:2 -> rejected by omitting -M
        # ambig1: R1 at pos 1300, R2 at pos 1400 (overlaps GeneB and GeneC) -> Unassigned_Ambiguity (omitting -O)
        sam_content = (
            "@HD\tVN:1.6\tSO:unsorted\n"
            "@SQ\tSN:chr1\tLN:10000\n"
            "@SQ\tSN:chr2\tLN:10000\n"
            # pair1 (valid)
            "pair1\t99\tchr1\t150\t60\t50M\t=\t300\t200\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            "pair1\t147\tchr1\t300\t60\t50M\t=\t150\t-200\t" + "T"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            # pair2 (valid)
            "pair2\t99\tchr1\t200\t60\t50M\t=\t350\t200\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            "pair2\t147\tchr1\t350\t60\t50M\t=\t200\t-200\t" + "T"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            # orphan1 (mate unmapped)
            "orphan1\t73\tchr1\t250\t60\t50M\t*\t0\t0\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            "orphan1\t133\t*\t0\t0\t*\tchr1\t250\t0\t" + "T"*50 + "\t" + "F"*50 + "\n"
            # chimeric1 (ends on different chromosomes)
            "chimeric1\t97\tchr1\t220\t60\t50M\tchr2\t200\t0\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            "chimeric1\t145\tchr2\t200\t60\t50M\tchr1\t220\t0\t" + "T"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            # multimap1 (NH:i:2)
            "multimap1\t99\tchr1\t260\t0\t50M\t=\t400\t190\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:2\n"
            "multimap1\t147\tchr1\t400\t0\t50M\t=\t260\t-190\t" + "T"*50 + "\t" + "F"*50 + "\tNH:i:2\n"
            # ambig1 (overlaps GeneB and GeneC)
            "ambig1\t99\tchr1\t1300\t60\t50M\t=\t1400\t150\t" + "A"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
            "ambig1\t147\tchr1\t1400\t60\t50M\t=\t1300\t-150\t" + "T"*50 + "\t" + "F"*50 + "\tNH:i:1\n"
        )
        sam_path.write_text(sam_content, encoding="utf-8")

        # Convert SAM to sorted BAM using WSL samtools
        # Convert paths for execution
        wsl_sam = to_wsl_path(sam_path)
        wsl_bam = to_wsl_path(bam_path)
        wsl_sorted = to_wsl_path(sorted_bam_path)
        wsl_gtf = to_wsl_path(gtf_path)
        wsl_fc = to_wsl_path(fc_out_path)

        res_sam = run_wsl(f"samtools view -bS {wsl_sam} > {wsl_bam} && samtools sort {wsl_bam} -o {wsl_sorted}")
        assert res_sam.returncode == 0, f"samtools failed: {res_sam.stderr}"

        # 3. Execute featureCounts with governed flags: -p --countReadPairs -B -C -t exon -g gene_id
        # (omitting -M and -O)
        fc_cmd = f"featureCounts -p --countReadPairs -B -C -t exon -g gene_id -s 0 -a {wsl_gtf} -o {wsl_fc} {wsl_sorted}"
        res_fc = run_wsl(fc_cmd)
        assert res_fc.returncode == 0, f"featureCounts failed: {res_fc.stderr}"

        # 4. Parse featureCounts results
        assert fc_out_path.exists(), "featureCounts output missing"
        counts = {}
        with fc_out_path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                parts = line.strip().split("\t")
                counts[parts[0]] = int(parts[6])

        # Assertions for Test A:
        # Expected: GeneA (ENSG00000000001) has count == 2 (pair1, pair2)
        # Fails if 4 (mate-level counting), fails if orphan or chimeric counted!
        assert counts.get("ENSG00000000001") == 2, (
            f"Test A failed: expected count 2 for GeneA, got {counts.get('ENSG00000000001')}. "
            "Accidental mate-level counting, orphan counting, or chimeric counting occurred!"
        )

        # Assertions for Test B:
        # Multimapper multimap1 has NH:i:2 and was not counted in GeneA
        # (Verified because count is exactly 2, not 3)

        # Assertions for Test C:
        # Ambiguous fragment ambig1 overlaps GeneB and GeneC.
        # With -O omitted, both GeneB and GeneC should have count == 0!
        assert counts.get("ENSG00000000002") == 0, f"Test C failed: GeneB has non-zero count {counts.get('ENSG00000000002')}"
        assert counts.get("ENSG00000000003") == 0, f"Test C failed: GeneC has non-zero count {counts.get('ENSG00000000003')}"

        # Check featureCounts summary file for Unassigned_Ambiguity and Unassigned_Unmapped
        summary_path = Path(str(fc_out_path) + ".summary")
        assert summary_path.exists()
        summary_text = summary_path.read_text(encoding="utf-8")
        assert "Unassigned_Ambiguity\t1" in summary_text, f"Expected 1 Unassigned_Ambiguity in summary:\n{summary_text}"
        assert "Assigned\t2" in summary_text, f"Expected exactly 2 Assigned fragments in summary:\n{summary_text}"


def test_d_adapter_numerical_conservation():
    """Execute Test D: Adapter numerical-conservation test with 0 discrepancies."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        raw_fc_file = tmp_p / "raw_featurecounts.txt"
        out_matrix = tmp_p / "prjna802072_pilot_raw_counts.tsv"
        out_universe = tmp_p / "prjna802072_featurecounts_gene_universe.tsv"

        # Construct synthetic featureCounts native output with 12 samples
        header_bams = "\t".join([f"/path/to/alignments/{run}_Aligned.sortedByCoord.out.bam" for run in APPROVED_RUN_ACCESSIONS])
        lines = [
            "# Program:featureCounts v2.0.6; Command:featureCounts ...",
            f"Geneid\tChr\tStart\tEnd\tStrand\tLength\t{header_bams}",
        ]

        expected_counts = {}
        for i in range(100):
            gid = f"ENSG{i:011d}"
            chr_val = f"chr{(i%22)+1}"
            start_val = f"{i*1000 + 1}"
            end_val = f"{i*1000 + 500}"
            strand_val = "+" if i % 2 == 0 else "-"
            len_val = "500"
            sample_vals = [(i * 13 + j * 7) for j in range(12)]
            expected_counts[gid] = sample_vals
            lines.append(f"{gid}\t{chr_val}\t{start_val}\t{end_val}\t{strand_val}\t{len_val}\t" + "\t".join(str(v) for v in sample_vals))

        raw_fc_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

        # Run adapter
        res = convert_featurecounts_to_matrix(
            raw_featurecounts_path=raw_fc_file,
            output_matrix_path=out_matrix,
            output_universe_path=out_universe,
            approved_runs=APPROVED_RUN_ACCESSIONS,
        )

        assert res["status"] == "CONVERSION_SUCCESS"
        assert res["numerical_conservation_discrepancies"] == 0
        assert res["is_numerically_conserved"] is True
        assert res["total_genes_emitted"] == 100

        # Verify output matrix content
        matrix_lines = out_matrix.read_text(encoding="utf-8").strip().split("\n")
        assert len(matrix_lines) == 101  # header + 100 genes
        header = matrix_lines[0].split("\t")
        assert header == ["gene_id"] + list(APPROVED_RUN_ACCESSIONS)

        for line in matrix_lines[1:]:
            parts = line.split("\t")
            gid = parts[0]
            counts = [int(x) for x in parts[1:]]
            assert counts == expected_counts[gid], f"Discrepancy for gene {gid}!"


def test_e_sample_header_mapping():
    """Execute Test E: Verify deterministic mapping of BAM column headers to approved Run accessions."""
    for run in APPROVED_RUN_ACCESSIONS:
        valid_headers = [
            f"alignments/{run}_Aligned.sortedByCoord.out.bam",
            f"/scratch/data/alignments/{run}.bam",
            f"{run}.bam",
            f"STAR_output/{run}_Aligned.out.bam",
        ]
        for vh in valid_headers:
            mapped = map_bam_header_to_run_id(vh, APPROVED_RUN_ACCESSIONS)
            assert mapped == run, f"Failed mapping '{vh}' -> expected '{run}', got '{mapped}'"

    # Test ambiguous header raises ValueError
    with pytest.raises(ValueError, match="Ambiguous"):
        map_bam_header_to_run_id(f"alignments/{APPROVED_RUN_ACCESSIONS[0]}_{APPROVED_RUN_ACCESSIONS[1]}.bam", APPROVED_RUN_ACCESSIONS)

    # Test unknown header raises ValueError
    with pytest.raises(ValueError, match="does not match"):
        map_bam_header_to_run_id("alignments/UNKNOWN_RUN_Aligned.bam", APPROVED_RUN_ACCESSIONS)


def test_f_bed12_reproducibility():
    """Execute Test F: Verify deterministic BED12 reference generation (identical SHA-256 and content)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        gtf_path = tmp_p / "test.gtf"
        bed_out_1 = tmp_p / "ref_run1.bed"
        bed_out_2 = tmp_p / "ref_run2.bed"

        gtf_content = (
            'chr1\tensembl\tgene\t100\t2000\t.\t+\t.\tgene_id "ENSG1";\n'
            'chr1\tensembl\ttranscript\t100\t2000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
            'chr1\tensembl\texon\t100\t300\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
            'chr1\tensembl\texon\t500\t800\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
            'chr1\tensembl\texon\t1200\t2000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
            'chr2\tensembl\tgene\t300\t1500\t.\t-\t.\tgene_id "ENSG2";\n'
            'chr2\tensembl\ttranscript\t300\t1500\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2";\n'
            'chr2\tensembl\texon\t300\t600\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2";\n'
            'chr2\tensembl\texon\t1000\t1500\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2";\n'
        )
        gtf_path.write_text(gtf_content, encoding="utf-8")

        res1 = generate_bed12_from_gtf(gtf_path, bed_out_1)
        res2 = generate_bed12_from_gtf(gtf_path, bed_out_2)

        assert res1["status"] == "BED12_GENERATION_SUCCESS"
        assert res2["status"] == "BED12_GENERATION_SUCCESS"
        assert res1["output_bed_sha256"] == res2["output_bed_sha256"], "Non-deterministic SHA-256 in BED12 output!"
        assert res1["total_transcripts_written"] == 2
        assert res2["total_transcripts_written"] == 2

        content1 = bed_out_1.read_text(encoding="utf-8")
        content2 = bed_out_2.read_text(encoding="utf-8")
        assert content1 == content2, "Content mismatch between successive runs!"


def test_g_frozen_pipeline_immutability():
    """Execute Test G: Verify frozen pipeline files have not been modified."""
    repo_root = Path(__file__).resolve().parents[2]
    
    # Check git diff for protected frozen paths
    protected_paths = [
        "src/qc",
        "src/gene_harmonization",
        "src/normalization",
        "src/differential_expression",
        "src/response_builder",
        "src/features/cross_study",
        "config/qc_rules.yaml",
        "config/gene_harmonization.yaml",
        "config/bulk_normalization.yaml",
        "data/references/ensembl_human_genes_v112.tsv",
    ]

    for p in protected_paths:
        res = subprocess.run(
            ["git", "diff", "--ignore-space-at-eol", "--exit-code", "HEAD", "--", p],
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0, f"Frozen path '{p}' has uncommitted modifications:\n{res.stdout}"
