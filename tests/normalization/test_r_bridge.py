"""Unit tests for the official DESeq2 R execution bridge."""

import json
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch
import pytest

from src.normalization.config import NormalizationConfig
from src.normalization.models import Severity
from src.normalization.r_bridge import (
    REnvironmentInfo,
    estimate_size_factors_r,
    get_environment_identity,
    probe_r_environment,
)


def test_probe_r_environment_when_binary_missing() -> None:
    info = probe_r_environment(r_binary="nonexistent_r_binary_12345")
    assert not info.is_available
    assert info.r_version == "missing"
    assert info.deseq2_version == "missing"
    assert "was not found" in (info.error_message or "")


def test_estimate_size_factors_r_fails_cleanly_when_r_missing() -> None:
    cfg = NormalizationConfig(r_binary_path="nonexistent_r_binary_12345")
    cols = ((10, 20), (30, 40))
    genes = ("ENSG001", "ENSG002")
    samples = ("s1", "s2")

    sfs, env_info, findings, prov = estimate_size_factors_r(
        matrix_columns=cols,
        gene_ids=genes,
        sample_ids=samples,
        config=cfg,
    )
    assert sfs is None
    assert not env_info.is_available
    assert any(f.rule_id == "r_environment_missing" for f in findings)
    assert any(f.severity == Severity.ERROR for f in findings)


def test_environment_identity_tracks_renv_lock() -> None:
    ident = get_environment_identity()
    assert ident.startswith("renv:") or ident == "unlocked_environment"


@patch("src.normalization.r_bridge.probe_r_environment")
@patch("subprocess.run")
def test_estimate_size_factors_r_success(mock_subproc: MagicMock, mock_probe: MagicMock, tmp_path: Path) -> None:
    mock_probe.return_value = REnvironmentInfo(
        r_version="4.4.3",
        bioc_version="3.20",
        deseq2_version="1.46.0",
        operating_system="Windows",
        execution_timestamp="2026-09-24T12:00:00Z",
        environment_identity="renv:test",
        is_available=True,
    )

    # Mock subprocess creating the output JSON
    def side_effect(cmd: list[str], **kwargs: Any) -> MagicMock:
        output_idx = cmd.index("--output")
        out_file = Path(cmd[output_idx + 1])
        out_file.write_text(
            json.dumps({
                "status": "SUCCESS",
                "method": "ratio",
                "size_factors": {"s1": 0.85, "s2": 1.15},
                "r_version": "4.4.3",
                "bioc_version": "3.20",
                "deseq2_version": "1.46.0",
                "operating_system": "Windows",
                "execution_timestamp": "2026-09-24T12:00:00Z",
                "error_message": None,
            })
        )
        return MagicMock(returncode=0, stdout="", stderr="")

    mock_subproc.side_effect = side_effect

    cols = ((10, 20), (30, 40))
    genes = ("ENSG001", "ENSG002")
    samples = ("s1", "s2")

    sfs, env_info, findings, prov = estimate_size_factors_r(
        matrix_columns=cols,
        gene_ids=genes,
        sample_ids=samples,
        r_script_path=Path("scripts/r/estimate_size_factors.R"),
    )
    assert sfs == {"s1": 0.85, "s2": 1.15}
    assert not any(f.severity == Severity.ERROR for f in findings)
    assert "input_matrix_sha256" in prov


@patch("src.normalization.r_bridge.probe_r_environment")
@patch("subprocess.run")
def test_estimate_size_factors_r_standard_zero_geomean_failure(mock_subproc: MagicMock, mock_probe: MagicMock) -> None:
    mock_probe.return_value = REnvironmentInfo(
        r_version="4.4.3",
        bioc_version="3.20",
        deseq2_version="1.46.0",
        operating_system="Windows",
        execution_timestamp="2026-09-24T12:00:00Z",
        environment_identity="renv:test",
        is_available=True,
    )

    def side_effect(cmd: list[str], **kwargs: Any) -> MagicMock:
        output_idx = cmd.index("--output")
        out_file = Path(cmd[output_idx + 1])
        out_file.write_text(
            json.dumps({
                "status": "ERROR",
                "method": "ratio",
                "size_factors": None,
                "r_version": "4.4.3",
                "bioc_version": "3.20",
                "deseq2_version": "1.46.0",
                "operating_system": "Windows",
                "execution_timestamp": "2026-09-24T12:00:00Z",
                "error_message": "every gene contains at least one zero, cannot compute log geometric means",
            })
        )
        return MagicMock(returncode=4, stdout="", stderr="")

    mock_subproc.side_effect = side_effect

    cols = ((0, 20), (30, 0))
    genes = ("ENSG001", "ENSG002")
    samples = ("s1", "s2")

    sfs, env_info, findings, prov = estimate_size_factors_r(
        matrix_columns=cols,
        gene_ids=genes,
        sample_ids=samples,
        r_script_path=Path("scripts/r/estimate_size_factors.R"),
    )
    assert sfs is None
    assert any(f.rule_id == "standard_size_factors_failed" for f in findings)
    assert any(f.severity == Severity.REVIEW for f in findings)


@patch("src.normalization.r_bridge.probe_r_environment")
def test_version_mismatch_surfacing(mock_probe: MagicMock) -> None:
    mock_probe.return_value = REnvironmentInfo(
        r_version="4.3.1",
        bioc_version="3.18",
        deseq2_version="1.42.0",
        operating_system="Windows",
        execution_timestamp="2026-09-24T12:00:00Z",
        environment_identity="unlocked",
        is_available=True,
    )

    cfg_warn = NormalizationConfig(strict_version_check=False)
    cols = ((10,), (20,))
    genes = ("ENSG001",)
    samples = ("s1", "s2")

    with patch("subprocess.run") as mock_subproc:
        def side_effect(cmd: list[str], **kwargs: Any) -> MagicMock:
            out_idx = cmd.index("--output")
            Path(cmd[out_idx + 1]).write_text(json.dumps({"status": "SUCCESS", "size_factors": {"s1": 1.0, "s2": 1.0}}))
            return MagicMock(returncode=0)
        mock_subproc.side_effect = side_effect

        sfs, _, findings, _ = estimate_size_factors_r(cols, genes, samples, config=cfg_warn)
        assert any(f.rule_id == "r_version_mismatch" and f.severity == Severity.WARNING for f in findings)

    cfg_strict = NormalizationConfig(strict_version_check=True)
    sfs, _, findings, _ = estimate_size_factors_r(cols, genes, samples, config=cfg_strict)
    assert sfs is None
    assert any(f.rule_id == "r_version_mismatch" and f.severity == Severity.ERROR for f in findings)
