"""Configuration loader, schema validator, and checksum calculation for Normalization v1."""

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any
import yaml

DEFAULT_CONFIG_PATH = Path("config/bulk_normalization.yaml")


@dataclass(frozen=True, slots=True)
class NormalizationConfig:
    """Configurable parameters for Bulk RNA-seq Normalization v1.

    All scientific thresholds default to None (unset) to eliminate invented thresholds.
    """

    size_factor_method: str = "ratio"
    min_evaluable_genes: int | None = None
    max_size_factor_ratio: float | None = None
    min_size_factor: float | None = None
    max_size_factor: float | None = None
    r_binary_path: str = "Rscript"
    r_timeout_seconds: int = 300
    expected_r_version: str = "4.4.3"
    expected_bioc_version: str = "3.20"
    expected_deseq2_version: str = "1.46.0"
    strict_version_check: bool = True
    config_version: str = "1.0.0"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "NormalizationConfig":
        valid_keys = {
            "size_factor_method",
            "min_evaluable_genes",
            "max_size_factor_ratio",
            "min_size_factor",
            "max_size_factor",
            "r_binary_path",
            "r_timeout_seconds",
            "expected_r_version",
            "expected_bioc_version",
            "expected_deseq2_version",
            "strict_version_check",
            "config_version",
        }
        filtered = {k: v for k, v in data.items() if k in valid_keys}
        return cls(**filtered)


def validate_config(config_dict: dict[str, Any]) -> None:
    """Validate normalization configuration semantics."""
    if not isinstance(config_dict, dict):
        raise ValueError("Configuration must be a dictionary.")

    method = config_dict.get("size_factor_method", "ratio")
    if method not in ("ratio", "poscounts"):
        raise ValueError(f"Unsupported size_factor_method: '{method}'. Must be 'ratio' or 'poscounts'.")

    for key in ("min_evaluable_genes", "max_size_factor_ratio", "min_size_factor", "max_size_factor"):
        val = config_dict.get(key)
        if val is not None:
            if not isinstance(val, (int, float)) or val < 0:
                raise ValueError(f"Configuration parameter '{key}' must be a non-negative number or None.")


def config_checksum(config: NormalizationConfig | dict[str, Any]) -> str:
    """Compute deterministic SHA-256 checksum of normalized configuration."""
    raw = config.to_dict() if isinstance(config, NormalizationConfig) else dict(config)
    serialized = json.dumps(raw, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def load_config(path: str | Path | None = None) -> NormalizationConfig:
    """Load NormalizationConfig from YAML path or return defaults if absent."""
    target_path = Path(path) if path else DEFAULT_CONFIG_PATH
    if not target_path.exists():
        return NormalizationConfig()

    with target_path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    validate_config(data)
    return NormalizationConfig.from_dict(data)
