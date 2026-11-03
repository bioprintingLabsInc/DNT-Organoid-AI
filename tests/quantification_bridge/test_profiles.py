"""Unit and preflight tests for execution profile governance and resource validation.

Verifies:
1. Portable execution profile parsing from YAML config.
2. Hardware baseline and tool inventory probing.
3. Preflight compatibility evaluation (fatal errors vs. advisory warnings).
4. Portable path resolution without machine-specific hardcoding.
5. QC metric threshold evaluation policies.
6. QuantificationBridgePipeline profile integration.
"""

import os
from pathlib import Path
import tempfile
import pytest

from src.quantification_bridge.profiles import (
    ExecutionProfile,
    HardwareBaseline,
    PreflightCheckResult,
    PreflightFinding,
    ToolInventory,
    evaluate_qc_metric,
    get_execution_profile,
    load_execution_profiles,
    probe_hardware_baseline,
    probe_tool_inventory,
    resolve_pipeline_paths,
    run_profile_preflight_checks,
    validate_pilot_configuration,
)
from src.quantification_bridge.pipeline import QuantificationBridgePipeline


def test_load_execution_profiles_schema():
    """Verify execution_profiles.yaml loads and contains all required profiles and sections."""
    cfg = load_execution_profiles()
    assert "profiles" in cfg
    assert "qc_thresholds" in cfg
    assert "sparse_star_benchmarking_plan" in cfg
    assert "software_inventory" in cfg
    assert "reference_workflow" in cfg

    profiles = cfg["profiles"]
    expected_profiles = {"institutional_hpc", "workstation_standard", "workstation_evaluation"}
    assert expected_profiles.issubset(set(profiles.keys()))

    for p_name in expected_profiles:
        p = profiles[p_name]
        assert "hardware_minimum" in p
        assert "star" in p
        assert "samtools" in p
        assert "reproducibility" in p
        assert "physical_ram_gb" in p["hardware_minimum"]
        assert "genome_sa_sparse_d" in p["star"]


def test_get_execution_profile_explicit_and_env():
    """Verify profile instantiation with explicit names and environment variable override."""
    # 1. Institutional HPC
    hpc = get_execution_profile("institutional_hpc")
    assert hpc.name == "institutional_hpc"
    assert hpc.physical_ram_gb_min >= 60.0
    assert hpc.genome_sa_sparse_d == 1
    assert hpc.external_sort_enabled is False

    # 2. Workstation Evaluation
    eval_prof = get_execution_profile("workstation_evaluation")
    assert eval_prof.name == "workstation_evaluation"
    assert eval_prof.physical_ram_gb_min <= 16.0
    assert eval_prof.genome_sa_sparse_d == 4
    assert eval_prof.genome_sa_index_nbases == 13
    assert eval_prof.external_sort_enabled is True
    assert eval_prof.samtools_sort_mem_per_thread == "1G"

    # 3. Invalid profile raises ValueError
    with pytest.raises(ValueError, match="Unknown execution profile"):
        get_execution_profile("non_existent_cluster")

    # 4. Environment variable override
    orig_env = os.environ.get("DNT_EXECUTION_PROFILE")
    try:
        os.environ["DNT_EXECUTION_PROFILE"] = "workstation_standard"
        env_prof = get_execution_profile()
        assert env_prof.name == "workstation_standard"
        assert env_prof.genome_sa_sparse_d == 2
    finally:
        if orig_env is not None:
            os.environ["DNT_EXECUTION_PROFILE"] = orig_env
        else:
            os.environ.pop("DNT_EXECUTION_PROFILE", None)


def test_resolve_pipeline_paths_portable():
    """Verify directory resolution handles default relative paths and environment overrides."""
    repo_root = Path(__file__).resolve().parents[2]
    paths = resolve_pipeline_paths(repo_root)

    assert paths["repo_root"] == repo_root
    assert paths["data_dir"] == repo_root / "data"
    assert paths["ref_dir"] == repo_root / "references" / "ensembl_112"
    assert "Homo_sapiens.GRCh38.112.gtf.gz" in str(paths["gtf_path"])

    # Test environment overrides
    with tempfile.TemporaryDirectory() as tmp_data, tempfile.TemporaryDirectory() as tmp_ref:
        orig_data = os.environ.get("DNT_DATA_DIR")
        orig_ref = os.environ.get("DNT_REF_DIR")
        try:
            os.environ["DNT_DATA_DIR"] = tmp_data
            os.environ["DNT_REF_DIR"] = tmp_ref
            custom_paths = resolve_pipeline_paths(repo_root)

            assert custom_paths["data_dir"] == Path(tmp_data).resolve()
            assert custom_paths["ref_dir"] == Path(tmp_ref).resolve()
            assert custom_paths["raw_fastqs_dir"].is_relative_to(Path(tmp_data).resolve())
            assert custom_paths["star_index_dir"].is_relative_to(Path(tmp_ref).resolve())
        finally:
            if orig_data is not None:
                os.environ["DNT_DATA_DIR"] = orig_data
            else:
                os.environ.pop("DNT_DATA_DIR", None)

            if orig_ref is not None:
                os.environ["DNT_REF_DIR"] = orig_ref
            else:
                os.environ.pop("DNT_REF_DIR", None)


def test_hardware_baseline_probe():
    """Verify hardware baseline empirical inspection."""
    hw = probe_hardware_baseline()
    assert isinstance(hw, HardwareBaseline)
    assert hw.cpu_cores >= 1
    assert hw.physical_ram_gb > 0.0
    assert hw.data_dir_free_gb >= 0.0
    assert hw.ref_dir_free_gb >= 0.0

    d = hw.to_dict()
    assert "cpu_cores" in d
    assert "physical_ram_gb" in d
    assert "swap_gb" in d
    assert "os_kernel" in d


def test_tool_inventory_probe():
    """Verify bioinformatics software inventory auditing."""
    inv = probe_tool_inventory()
    assert isinstance(inv, ToolInventory)
    required = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py", "fasterq-dump"]
    for t in required:
        assert t in inv.tools
        t_info = inv.tools[t]
        assert "path" in t_info
        assert "version" in t_info


def test_run_profile_preflight_checks_ram_enforcement():
    """Verify that resource checks distinguish fatal blockers from advisory warnings."""
    eval_prof = get_execution_profile("workstation_evaluation")
    res_eval = run_profile_preflight_checks(eval_prof, enforce_physical_ram=False)
    assert isinstance(res_eval, PreflightCheckResult)

    hpc_prof = get_execution_profile("institutional_hpc")

    # If host has < 64 GB RAM, enforce_physical_ram=True must yield FATAL
    res_hpc_enforced = run_profile_preflight_checks(hpc_prof, enforce_physical_ram=True)
    if res_hpc_enforced.hardware.physical_ram_gb < 64.0:
        assert res_hpc_enforced.has_fatal_errors
        assert not res_hpc_enforced.is_compatible
        ram_findings = [f for f in res_hpc_enforced.findings if f.code == "insufficient_physical_ram"]
        assert len(ram_findings) >= 1
        assert ram_findings[0].severity == "FATAL"

    # With enforce_physical_ram=False, insufficient RAM is logged as WARNING, not FATAL
    res_hpc_advisory = run_profile_preflight_checks(hpc_prof, enforce_physical_ram=False)
    if res_hpc_advisory.hardware.physical_ram_gb < 64.0:
        ram_warnings = [f for f in res_hpc_advisory.findings if f.code == "ram_below_profile_minimum"]
        assert len(ram_warnings) >= 1
        assert ram_warnings[0].severity == "WARNING"


def test_evaluate_qc_metric_policies():
    """Verify multi-tier QC threshold evaluation: PASS, WARNING, FATAL."""
    # 1. Fatal blockers
    status, msg = evaluate_qc_metric("mapped_reads", 0.0)
    assert status == "FATAL"
    assert "Zero mapped reads" in msg

    status, msg = evaluate_qc_metric("fastq_file_size", 0.0)
    assert status == "FATAL"
    assert "Empty FASTQ" in msg

    # 2. Advisory warnings
    status, msg = evaluate_qc_metric("star_unique_mapping_fraction", 0.72)
    assert status == "WARNING"
    assert "unique mapping fraction" in msg

    status, msg = evaluate_qc_metric("star_multimapping_fraction", 0.22)
    assert status == "WARNING"
    assert "multimapping fraction" in msg

    status, msg = evaluate_qc_metric("fastp_q30_bases_fraction", 0.81)
    assert status == "WARNING"
    assert "Q30 base fraction" in msg

    status, msg = evaluate_qc_metric("featurecounts_assigned_fraction", 0.52)
    assert status == "WARNING"
    assert "assigned exonic read fraction" in msg

    # 3. Passing metrics
    status, msg = evaluate_qc_metric("star_unique_mapping_fraction", 0.88)
    assert status == "PASS"

    status, msg = evaluate_qc_metric("fastp_q30_bases_fraction", 0.94)
    assert status == "PASS"


def test_quantification_bridge_pipeline_profile_integration():
    """Verify QuantificationBridgePipeline initializes with designated profile and executes checks."""
    pipeline = QuantificationBridgePipeline(profile_name="workstation_evaluation")
    assert pipeline.profile.name == "workstation_evaluation"
    assert pipeline.profile.genome_sa_sparse_d == 4
    assert pipeline.profile.external_sort_enabled is True

    # Preflight check execution
    preflight = pipeline.run_profile_checks(enforce_physical_ram=False)
    assert isinstance(preflight, PreflightCheckResult)
    assert preflight.profile.name == "workstation_evaluation"


def test_validate_pilot_configuration_valid():
    """Assert validate_pilot_configuration succeeds with current configuration."""
    res = validate_pilot_configuration()
    assert res["status"] == "PILOT_CONFIGURATION_VALID"
    assert res["control_run"] == "SRR17834920"
    assert res["treatment_run"] == "SRR17834905"
    assert res["control_dose_um"] == 0.0
    assert res["treatment_dose_um"] == 250.0


def test_preflight_safeguards_memory_headroom_and_storage():
    """Assert preflight checks evaluate memory headroom and peak disk requirements."""
    eval_prof = get_execution_profile("workstation_evaluation")
    res = run_profile_preflight_checks(eval_prof, enforce_physical_ram=False)
    
    # Assert findings codes are structured properly
    finding_codes = {f.code for f in res.findings}
    # Should not have fatal errors on current setup if tools and space are present
    assert not any(f.severity == "FATAL" for f in res.findings if f.code == "insufficient_disk_space_pilot")
    
    # Provenance contains complete breakdown
    assert "profile" in res.provenance
    assert "hardware" in res.provenance
    assert "tools" in res.provenance

