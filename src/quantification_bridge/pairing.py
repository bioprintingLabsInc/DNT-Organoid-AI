"""Deterministic FASTQ pair grouping and validation contract.

Governs pairing of forward (R1) and reverse (R2) FASTQ reads by canonical sample identity,
preventing orphan reads, duplicate assignments, and filesystem ordering dependencies.
"""

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any

from src.ingestion.models import RegisteredAsset


class FastqPairingError(ValueError):
    """Raised when paired-end FASTQ grouping or integrity validation fails."""


@dataclass(frozen=True, slots=True)
class FastqPair:
    """Immutable record of a validated paired-end FASTQ library for a single sample."""

    sample_id: str
    r1_path: Path
    r2_path: Path
    assay_id: str | None = None
    r1_asset_id: str | None = None
    r2_asset_id: str | None = None
    extra_metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "r1_path": str(self.r1_path),
            "r2_path": str(self.r2_path),
            "assay_id": self.assay_id,
            "r1_asset_id": self.r1_asset_id,
            "r2_asset_id": self.r2_asset_id,
            "extra_metadata": self.extra_metadata,
        }


# Regular expressions for identifying read mates in filenames
# Matches patterns like:
# _1.fastq.gz, _2.fastq.gz, _1.fq.gz, _2.fq.gz
# _R1.fastq.gz, _R2.fastq.gz, _R1_001.fastq.gz, _R2_001.fastq.gz
# .1.fastq.gz, .2.fastq.gz
_R1_PATTERNS = [
    re.compile(r"[_.-]R?1(_\d+)?\.(fastq|fq)(\.gz)?$", re.IGNORECASE),
    re.compile(r"[_.-]1\.(fastq|fq)(\.gz)?$", re.IGNORECASE),
    re.compile(r"[_.-]R1(_\d+)?$", re.IGNORECASE),
]
_R2_PATTERNS = [
    re.compile(r"[_.-]R?2(_\d+)?\.(fastq|fq)(\.gz)?$", re.IGNORECASE),
    re.compile(r"[_.-]2\.(fastq|fq)(\.gz)?$", re.IGNORECASE),
    re.compile(r"[_.-]R2(_\d+)?$", re.IGNORECASE),
]


def determine_read_mate(filename_or_path: str | Path) -> str | None:
    """Determine whether a file is Read 1 (R1) or Read 2 (R2).
    
    Returns 'R1', 'R2', or None if mate cannot be unambiguously resolved.
    """
    name = Path(filename_or_path).name
    is_r1 = any(p.search(name) for p in _R1_PATTERNS)
    is_r2 = any(p.search(name) for p in _R2_PATTERNS)

    if is_r1 and not is_r2:
        return "R1"
    if is_r2 and not is_r1:
        return "R2"
    return None


def extract_sample_stem(filename_or_path: str | Path) -> str:
    """Extract sample identifier prefix from a FASTQ filename by stripping mate designations."""
    name = Path(filename_or_path).name
    for p in _R1_PATTERNS + _R2_PATTERNS:
        subbed = p.sub("", name)
        if subbed != name:
            return subbed
    # Fallback: remove extensions
    stem = name
    for ext in (".gz", ".fastq", ".fq"):
        if stem.lower().endswith(ext):
            stem = stem[:-len(ext)]
    return stem


def group_fastq_pairs(
    assets: tuple[RegisteredAsset, ...],
    metadata: dict[str, Any] | None = None,
) -> tuple[FastqPair, ...]:
    """Deterministically group registered FASTQ assets into validated FastqPair objects.

    Validation Rules:
    1. Every FASTQ asset must unambiguously identify as R1 or R2.
    2. Every sample must have exactly one R1 and exactly one R2.
    3. No orphan mates (R1 without R2, or R2 without R1) are permitted.
    4. No duplicate assignments (multiple R1s or multiple R2s for the same sample).
    5. No FASTQ file may be assigned to multiple samples.
    6. Output order is deterministic (sorted by canonical sample_id), independent of filesystem order.

    Sample resolution priority:
    - If metadata sequencing_assays are present, maps asset.assay_ids -> assay -> sample_id.
    - Otherwise, derives sample identity from the filename stem.
    """
    if not assets:
        raise FastqPairingError("Cannot group FASTQ pairs: asset list is empty.")

    # 1. Build assay -> sample_id mapping from metadata if available
    assay_to_sample: dict[str, str] = {}
    if metadata and isinstance(metadata.get("sequencing_assays"), list):
        for assay in metadata["sequencing_assays"]:
            if isinstance(assay, dict) and "assay_id" in assay:
                aid = assay["assay_id"]
                sids = assay.get("sample_ids", [])
                if isinstance(sids, list) and len(sids) == 1:
                    assay_to_sample[aid] = str(sids[0])

    # 2. Inspect each asset
    # Map: sample_id -> {"R1": (asset, path), "R2": (asset, path)}
    sample_buckets: dict[str, dict[str, tuple[RegisteredAsset, Path]]] = {}
    seen_paths: dict[Path, str] = {}

    for asset in assets:
        local_path = asset.local_path
        resolved_path = local_path.resolve() if local_path.exists() else local_path

        # Check for single file assigned to multiple assets
        if resolved_path in seen_paths:
            raise FastqPairingError(
                f"FASTQ file '{local_path.name}' is assigned multiple times (assets '{seen_paths[resolved_path]}' and '{asset.asset_id}')."
            )
        seen_paths[resolved_path] = asset.asset_id

        # Determine mate
        mate = determine_read_mate(local_path)
        if mate is None:
            raise FastqPairingError(
                f"Cannot determine read mate (R1 vs R2) for asset '{asset.asset_id}' with file '{local_path.name}'."
            )

        # Determine sample_id
        sample_id: str | None = None
        target_assay_id: str | None = None
        if asset.assay_ids:
            for aid in asset.assay_ids:
                if aid in assay_to_sample:
                    sample_id = assay_to_sample[aid]
                    target_assay_id = aid
                    break

        if not sample_id:
            sample_id = extract_sample_stem(local_path)

        if not sample_id:
            raise FastqPairingError(
                f"Cannot establish sample identity for asset '{asset.asset_id}' with file '{local_path.name}'."
            )

        bucket = sample_buckets.setdefault(sample_id, {})
        if mate in bucket:
            existing_asset = bucket[mate][0]
            raise FastqPairingError(
                f"Duplicate {mate} assignment for sample '{sample_id}': assets '{existing_asset.asset_id}' and '{asset.asset_id}'."
            )
        bucket[mate] = (asset, local_path)

    # 3. Validate complete pairs and check for orphans
    pairs: list[FastqPair] = []
    for sample_id in sorted(sample_buckets.keys()):
        bucket = sample_buckets[sample_id]
        has_r1 = "R1" in bucket
        has_r2 = "R2" in bucket

        if has_r1 and not has_r2:
            r1_name = bucket["R1"][1].name
            raise FastqPairingError(
                f"Orphan mate detected for sample '{sample_id}': R1 '{r1_name}' has no matching R2 reverse read."
            )
        if has_r2 and not has_r1:
            r2_name = bucket["R2"][1].name
            raise FastqPairingError(
                f"Orphan mate detected for sample '{sample_id}': R2 '{r2_name}' has no matching R1 forward read."
            )

        r1_asset, r1_path = bucket["R1"]
        r2_asset, r2_path = bucket["R2"]

        # Reconcile assay ID
        shared_assays = tuple(set(r1_asset.assay_ids).intersection(set(r2_asset.assay_ids)))
        assay_id = shared_assays[0] if shared_assays else (r1_asset.assay_ids[0] if r1_asset.assay_ids else None)

        pairs.append(
            FastqPair(
                sample_id=sample_id,
                r1_path=r1_path,
                r2_path=r2_path,
                assay_id=assay_id,
                r1_asset_id=r1_asset.asset_id,
                r2_asset_id=r2_asset.asset_id,
            )
        )

    return tuple(pairs)
