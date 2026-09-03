"""Bounded, non-destructive inspection of untrusted local assets."""

from pathlib import Path

from .models import AssetFinding, Inspection


MAX_INSPECTION_BYTES = 262_144


def inspect_asset(path: Path) -> Inspection:
    """Read at most MAX_INSPECTION_BYTES and return structured failures."""
    suffixes = tuple(suffix.lower() for suffix in path.suffixes)
    try:
        if not path.is_file():
            return Inspection(False, False, None, b"", suffixes, (
                AssetFinding("ERROR", "missing_file", "source_reference", f"File does not exist: {path}"),
            ))
        size = path.stat().st_size
        with path.open("rb") as stream:
            prefix = stream.read(MAX_INSPECTION_BYTES)
    except (OSError, PermissionError) as error:
        return Inspection(path.exists(), False, None, b"", suffixes, (
            AssetFinding("ERROR", "read_failure", "source_reference", f"File cannot be inspected: {error}"),
        ))
    findings = () if size else (
        AssetFinding("ERROR", "empty_file", "source_reference", "File is empty."),
    )
    return Inspection(True, True, size, prefix, suffixes, findings)

