"""Unit tests for Normalization v1 configuration module."""

from pathlib import Path
import pytest
import yaml

from src.normalization.config import (
    NormalizationConfig,
    config_checksum,
    load_config,
    validate_config,
)


def test_default_config_thresholds_are_none() -> None:
    cfg = NormalizationConfig()
    assert cfg.size_factor_method == "ratio"
    assert cfg.min_evaluable_genes is None
    assert cfg.max_size_factor_ratio is None
    assert cfg.min_size_factor is None
    assert cfg.max_size_factor is None
    assert cfg.expected_r_version == "4.4.3"
    assert cfg.expected_bioc_version == "3.20"
    assert cfg.expected_deseq2_version == "1.46.0"


def test_validate_config_rejects_invalid_method() -> None:
    with pytest.raises(ValueError, match="Unsupported size_factor_method"):
        validate_config({"size_factor_method": "invalid_method"})


def test_validate_config_rejects_negative_thresholds() -> None:
    with pytest.raises(ValueError, match="must be a non-negative number"):
        validate_config({"min_evaluable_genes": -5})


def test_config_checksum_deterministic() -> None:
    cfg1 = NormalizationConfig()
    cfg2 = NormalizationConfig()
    assert config_checksum(cfg1) == config_checksum(cfg2)
    assert len(config_checksum(cfg1)) == 64


def test_load_config_from_yaml(tmp_path: Path) -> None:
    cfg_file = tmp_path / "custom_config.yaml"
    cfg_file.write_text(
        yaml.dump({
            "size_factor_method": "poscounts",
            "min_evaluable_genes": 50,
            "max_size_factor_ratio": 15.0,
        })
    )
    loaded = load_config(cfg_file)
    assert loaded.size_factor_method == "poscounts"
    assert loaded.min_evaluable_genes == 50
    assert loaded.max_size_factor_ratio == 15.0
