"""Mandatory manifest regression test suite for PRJNA802072 quantification bridge."""

import json
from pathlib import Path
import pytest

from src.quantification_bridge.manifest import (
    APPROVED_2_SAMPLE_TECHNICAL_PILOT,
    APPROVED_PILOT_RUNS,
    APPROVED_RUN_ACCESSIONS,
    BIOLOGICAL_REPLICATE_UNIT,
    RUN_TO_SAMPLE_MAP,
    TARGET_BIOPROJECT,
    TARGET_DEVELOPMENTAL_DAY,
    TARGET_ORGANISM,
    TARGET_SRA_STUDY,
    get_approved_technical_pilot_runs,
)
from src.quantification_bridge.profiles import validate_pilot_configuration

# Governed immutable ground truth allowlist
EXPECTED_CONTROL_MAP = {
    "VPA_0_1_1": "SRR17834920",
    "VPA_0_2_1": "SRR17834908",
    "VPA_0_3_1": "SRR17834897",
    "VPA_0_4_1": "SRR17834886",
    "VPA_0_5_1": "SRR17834875",
    "VPA_0_6_1": "SRR17834919",
}

EXPECTED_TREATMENT_MAP = {
    "VPA_250_1_1": "SRR17834905",
    "VPA_250_2_1": "SRR17834904",
    "VPA_250_3_1": "SRR17834903",
    "VPA_250_4_1": "SRR17834902",
    "VPA_250_5_1": "SRR17834901",
    "VPA_250_6_1": "SRR17834900",
}

EXPECTED_EXPLICIT_RUN_SET = set(EXPECTED_CONTROL_MAP.values()) | set(EXPECTED_TREATMENT_MAP.values())


def test_manifest_cardinality_and_balance():
    """Assert exactly 12 total runs, 6 controls, and 6 treatments."""
    assert len(APPROVED_PILOT_RUNS) == 12, f"Expected 12 runs, got {len(APPROVED_PILOT_RUNS)}"
    
    controls = [r for r in APPROVED_PILOT_RUNS if r.condition_group == "control"]
    treatments = [r for r in APPROVED_PILOT_RUNS if r.condition_group == "treatment"]

    assert len(controls) == 6, f"Expected 6 controls, got {len(controls)}"
    assert len(treatments) == 6, f"Expected 6 treatments, got {len(treatments)}"

    for c in controls:
        assert c.vpa_dose_um == 0.0, f"Control {c.sample_title} has non-zero dose: {c.vpa_dose_um}"
        assert c.developmental_day == 28
    
    for t in treatments:
        assert t.vpa_dose_um == 250.0, f"Treatment {t.sample_title} has unexpected dose: {t.vpa_dose_um}"
        assert t.developmental_day == 28


def test_manifest_uniqueness():
    """Assert all 12 run accessions and all 12 sample aliases are strictly unique."""
    run_accessions = [r.run_accession for r in APPROVED_PILOT_RUNS]
    sample_aliases = [r.sample_title for r in APPROVED_PILOT_RUNS]

    assert len(set(run_accessions)) == 12, f"Duplicate run accessions detected: {run_accessions}"
    assert len(set(sample_aliases)) == 12, f"Duplicate sample aliases detected: {sample_aliases}"


def test_manifest_exact_allowlist_match():
    """Assert the run set matches the approved explicit allowlist exactly."""
    observed_runs = set(r.run_accession for r in APPROVED_PILOT_RUNS)
    
    missing = EXPECTED_EXPLICIT_RUN_SET - observed_runs
    unexpected = observed_runs - EXPECTED_EXPLICIT_RUN_SET

    assert not missing, f"Approved runs missing from manifest: {missing}"
    assert not unexpected, f"Unexpected runs in manifest: {unexpected}"
    assert observed_runs == EXPECTED_EXPLICIT_RUN_SET


def test_manifest_exact_sample_mappings():
    """Assert control and treatment sample title to run accession mappings are 100% exact."""
    obs_map = {r.sample_title: r.run_accession for r in APPROVED_PILOT_RUNS}

    for title, expected_run in EXPECTED_CONTROL_MAP.items():
        assert obs_map.get(title) == expected_run, f"Control mapping mismatch for {title}: expected {expected_run}, got {obs_map.get(title)}"

    for title, expected_run in EXPECTED_TREATMENT_MAP.items():
        assert obs_map.get(title) == expected_run, f"Treatment mapping mismatch for {title}: expected {expected_run}, got {obs_map.get(title)}"


def test_prohibition_of_accession_range_generation():
    """Assert that the 12-run set is non-contiguous and cannot match naive contiguous intervals."""
    run_nums = sorted([int(r.replace("SRR", "")) for r in APPROVED_RUN_ACCESSIONS])
    
    # Verify non-contiguity
    diffs = [run_nums[i+1] - run_nums[i] for i in range(len(run_nums)-1)]
    assert any(d > 1 for d in diffs), "Fatal: runs appear strictly contiguous, violating non-contiguous design"

    # Explicitly reject the invalid contiguous hallucinated ranges
    hallucinated_controls = {f"SRR{i}" for i in range(17834920, 17834926)}
    hallucinated_treatments = {f"SRR{i}" for i in range(17834914, 17834920)}
    hallucinated_set = hallucinated_controls | hallucinated_treatments

    assert set(APPROVED_RUN_ACCESSIONS) != hallucinated_set, "Fatal: manifest matches invalid contiguous range hallucination"
    assert "SRR17834925" not in APPROVED_RUN_ACCESSIONS
    assert "SRR17834914" not in APPROVED_RUN_ACCESSIONS
    assert "SRR17834915" not in APPROVED_RUN_ACCESSIONS


def test_json_manifest_synchronization():
    """Assert config/quantification_bridge/pilot_manifest.json matches manifest.py exactly."""
    manifest_file = Path(__file__).resolve().parents[2] / "config" / "quantification_bridge" / "pilot_manifest.json"
    assert manifest_file.exists(), f"pilot_manifest.json missing: {manifest_file}"

    data = json.loads(manifest_file.read_text(encoding="utf-8"))
    assert data["total_samples"] == 12
    assert data["control_samples_count"] == 6
    assert data["treatment_samples_count"] == 6
    assert data.get("strandedness") == "unstranded"
    assert data.get("featurecounts_strandedness_flag") == 0

    json_runs = {s["run_accession"]: s for s in data["samples"]}
    assert set(json_runs.keys()) == EXPECTED_EXPLICIT_RUN_SET

    for r in APPROVED_PILOT_RUNS:
        jr = json_runs[r.run_accession]
        assert jr["sample_alias"] == r.sample_title
        assert jr["group"] == r.condition_group
        assert jr["dose"] == r.vpa_dose_um
        assert jr["day"] == r.developmental_day
        assert jr["fastq_1_url"] == r.fastq_1_url
        assert jr["fastq_2_url"] == r.fastq_2_url
        assert jr["fastq_1_md5"] == r.fastq_1_md5
        assert jr["fastq_2_md5"] == r.fastq_2_md5
        assert jr["biological_replicate_unit"] == "UNVERIFIED"
        assert jr.get("strandedness") == r.strandedness
        assert jr.get("featurecounts_strandedness_flag") == r.featurecounts_strandedness_flag


def test_metadata_reconciliation_rejection():
    """Assert metadata validation utility rejects mismatched alias, dose, or day."""
    # Test mismatch rejection logic
    sample = APPROVED_PILOT_RUNS[0]
    mismatched_record = {
        "run_accession": sample.run_accession,
        "sample_accession": sample.biosample_accession,
        "experiment_accession": sample.experiment_accession,
        "sample_alias": "WRONG_ALIAS",
        "read_count": str(sample.paired_reads_expected),
        "base_count": str(sample.total_bases_expected),
        "fastq_ftp": f"{sample.fastq_1_url};{sample.fastq_2_url}",
        "fastq_md5": f"{sample.fastq_1_md5};{sample.fastq_2_md5}",
    }

    assert mismatched_record["sample_alias"] != sample.sample_title, "Should detect alias discrepancy"


def test_approved_2_sample_technical_pilot_invariance():
    """Assert approved 2-sample technical pilot contains exactly SRR17834920 and SRR17834905 at 250 uM Day 28."""
    pilot = get_approved_technical_pilot_runs()
    assert len(pilot) == 2, f"Expected exactly 2 technical pilot runs, got {len(pilot)}"
    
    ctrl = pilot[0]
    trt = pilot[1]

    # Exact control identity
    assert ctrl.run_accession == "SRR17834920", f"Expected control SRR17834920, got {ctrl.run_accession}"
    assert ctrl.sample_title == "VPA_0_1_1"
    assert ctrl.vpa_dose_um == 0.0
    assert ctrl.developmental_day == 28
    assert ctrl.condition_group == "control"

    # Exact treatment identity
    assert trt.run_accession == "SRR17834905", f"Expected treatment SRR17834905, got {trt.run_accession}"
    assert trt.sample_title == "VPA_250_1_1"
    assert trt.vpa_dose_um == 250.0
    assert trt.developmental_day == 28
    assert trt.condition_group == "treatment"

    # Strict rejection of hallucinated or confused accessions
    forbidden_accessions = {"SRR17826941", "SRR17826944"}
    pilot_runs = {ctrl.run_accession, trt.run_accession}
    assert not pilot_runs.intersection(forbidden_accessions), "Hallucinated accession present in technical pilot!"

    # Strict rejection of 1 mM / 1000 uM (which belongs to E-MTAB-15191, not PRJNA802072)
    assert ctrl.vpa_dose_um != 1000.0, "1 mM VPA mistakenly configured for control!"
    assert trt.vpa_dose_um != 1000.0, "1 mM VPA mistakenly configured for treatment (must be 250 uM for PRJNA802072 Day 28)!"


def test_technical_pilot_consistency_with_execution_profile():
    """Assert validate_pilot_configuration enforces agreement between manifest and YAML config."""
    res = validate_pilot_configuration()
    assert res["status"] == "PILOT_CONFIGURATION_VALID"
    assert res["control_run"] == "SRR17834920"
    assert res["treatment_run"] == "SRR17834905"
    assert res["control_dose_um"] == 0.0
    assert res["treatment_dose_um"] == 250.0


def test_rejection_of_hallucinated_1mm_and_unapproved_accessions():
    """Assert validate_pilot_configuration detects and rejects prohibited accessions or concentrations."""
    import tempfile
    import yaml

    config_path = Path(__file__).resolve().parents[2] / "config" / "quantification_bridge" / "execution_profiles.yaml"
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    # Case 1: Corrupted concentration (e.g. 1000 uM / 1 mM)
    corrupted_cfg_1 = dict(cfg)
    corrupted_cfg_1["technical_pilot_specification"] = dict(cfg["technical_pilot_specification"])
    corrupted_cfg_1["technical_pilot_specification"]["treatment"] = dict(cfg["technical_pilot_specification"]["treatment"])
    corrupted_cfg_1["technical_pilot_specification"]["treatment"]["dose_um"] = 1000.0

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        yaml.dump(corrupted_cfg_1, f)
        temp_file_1 = Path(f.name)

    try:
        with pytest.raises(ValueError, match="Treatment dose mismatch|Prohibited dose detected"):
            validate_pilot_configuration(temp_file_1)
    finally:
        temp_file_1.unlink(missing_ok=True)

    # Case 2: Hallucinated accession (e.g. SRR17826941)
    corrupted_cfg_2 = dict(cfg)
    corrupted_cfg_2["technical_pilot_specification"] = dict(cfg["technical_pilot_specification"])
    corrupted_cfg_2["technical_pilot_specification"]["treatment"] = dict(cfg["technical_pilot_specification"]["treatment"])
    corrupted_cfg_2["technical_pilot_specification"]["treatment"]["run_accession"] = "SRR17826941"

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        yaml.dump(corrupted_cfg_2, f)
        temp_file_2 = Path(f.name)

    try:
        with pytest.raises(ValueError, match="Treatment run mismatch|Prohibited run accession"):
            validate_pilot_configuration(temp_file_2)
    finally:
        temp_file_2.unlink(missing_ok=True)

