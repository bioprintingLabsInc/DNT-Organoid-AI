"""Versioned QC rule loading and canonical checksums."""
import hashlib, json
from pathlib import Path
from typing import Any
import yaml

DEFAULT_RULES_PATH = Path(__file__).resolve().parents[2] / "config" / "qc_rules.yaml"
_TOP = {"rules_version", "ruleset_name", "scope", "structural_rules", "warning_thresholds", "outlier_rules", "replicate_correlation"}

def load_rules(path: str | Path = DEFAULT_RULES_PATH) -> dict[str, Any]:
    rules = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(rules, dict): raise ValueError("QC rules must be a mapping")
    unknown = set(rules) - _TOP
    if unknown: raise ValueError(f"Unknown QC rule keys: {sorted(unknown)}")
    for key in _TOP:
        if key not in rules: raise ValueError(f"Missing QC rules section: {key}")
    if rules["replicate_correlation"].get("method") != "spearman": raise ValueError("Replicate correlation method must be spearman")
    return rules

def rules_checksum(rules: dict[str, Any]) -> str:
    canonical = json.dumps(rules, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(canonical.encode()).hexdigest()
