"""Versioned configuration loading and canonical checksums for gene harmonization."""

import hashlib
import json
from pathlib import Path
from typing import Any
import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "gene_harmonization.yaml"

_REQUIRED_SECTIONS = {
    "config_version",
    "stage_name",
    "scope",
    "reference",
    "supported_identifier_types",
    "structural_rules",
    "coverage_thresholds",
}


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    """Validate gene harmonization configuration dictionary."""
    if not isinstance(config, dict):
        raise ValueError("Gene harmonization configuration must be a mapping.")

    missing = _REQUIRED_SECTIONS - set(config)
    if missing:
        raise ValueError(f"Missing required configuration sections: {sorted(missing)}")

    unknown = set(config) - _REQUIRED_SECTIONS
    if unknown:
        raise ValueError(f"Unknown configuration sections: {sorted(unknown)}")

    if config.get("stage_name") != "gene_harmonization":
        raise ValueError(f"Invalid stage_name '{config.get('stage_name')}'; expected 'gene_harmonization'.")

    return config


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load and validate the gene harmonization configuration."""
    target_path = Path(path) if path else DEFAULT_CONFIG_PATH
    if not target_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {target_path}")

    with target_path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)

    return validate_config(config)


def config_checksum(config: dict[str, Any]) -> str:
    """Compute a deterministic SHA-256 checksum of the configuration."""
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
