"""Structural validation for paired and single FASTQ sequencing files."""

import gzip
from pathlib import Path
from typing import Any


def _open_fastq(path: Path):
    """Open FASTQ file with transparent gzip decompression if file ends in .gz."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def validate_single_fastq(
    fastq_path: str | Path,
    sample_records_to_check: int = 100000,
) -> dict[str, Any]:
    """Validate a single-end FASTQ file for integrity and valid format.

    Checks:
    - File exists and is non-empty.
    - Successfully decompresses if gzipped.
    - Valid 4-line FASTQ structure (@header, seq, +, qual).
    - Sequence and quality strings have identical lengths.
    """
    p = Path(fastq_path)

    if not p.exists():
        raise FileNotFoundError(f"FASTQ file not found: {p}")

    size = p.stat().st_size
    if size == 0:
        raise ValueError(f"Empty FASTQ file detected: {p.name}")

    total_reads = 0
    checked_reads = 0

    try:
        with _open_fastq(p) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                p_plus = f.readline()
                q = f.readline()

                total_reads += 1

                if checked_reads < sample_records_to_check:
                    if not h.startswith("@"):
                        raise ValueError(f"Corrupt FASTQ header at record {total_reads}: '{h.strip()}'")
                    if not p_plus.startswith("+"):
                        raise ValueError(f"Corrupt FASTQ '+' delimiter at record {total_reads}")
                    if len(s.strip()) != len(q.strip()):
                        raise ValueError(f"Sequence and quality length mismatch at record {total_reads}")
                    checked_reads += 1

    except gzip.BadGzipFile as e:
        raise ValueError(f"Corrupted gzip stream in {p.name}: {e}")

    return {
        "status": "VALIDATION_PASSED",
        "fastq": str(p),
        "total_reads": total_reads,
        "records_structurally_verified": checked_reads,
        "read_layout": "single_end",
        "is_compressed": str(p).lower().endswith(".gz"),
    }


def validate_paired_fastqs(
    fastq_1_path: str | Path,
    fastq_2_path: str | Path,
    sample_records_to_check: int = 100000,
) -> dict[str, Any]:
    """Validate paired FASTQ files for synchronization, integrity, and valid format.

    Checks:
    - Files exist and are non-empty.
    - Successfully decompresses if gzipped.
    - Synchronized read headers between Read 1 and Read 2.
    - Valid 4-line FASTQ structure (@header, seq, +, qual).
    - Sequence and quality strings have identical lengths.
    - Total read pair counts agree.
    """
    p1 = Path(fastq_1_path)
    p2 = Path(fastq_2_path)

    if not p1.exists():
        raise FileNotFoundError(f"Read 1 FASTQ not found: {p1}")
    if not p2.exists():
        raise FileNotFoundError(f"Read 2 FASTQ not found: {p2}")

    size1 = p1.stat().st_size
    size2 = p2.stat().st_size
    if size1 == 0 or size2 == 0:
        raise ValueError(f"Empty FASTQ file detected: R1 size={size1}, R2 size={size2}")

    total_pairs = 0
    checked_pairs = 0

    try:
        with _open_fastq(p1) as f1, _open_fastq(p2) as f2:
            while True:
                h1 = f1.readline()
                h2 = f2.readline()

                if not h1 and not h2:
                    break
                if bool(h1) != bool(h2):
                    raise ValueError(f"Desynchronized EOF: R1 and R2 have unequal read counts in {p1.name} / {p2.name}")

                s1 = f1.readline()
                s2 = f2.readline()
                p_plus1 = f1.readline()
                p_plus2 = f2.readline()
                q1 = f1.readline()
                q2 = f2.readline()

                total_pairs += 1

                # Structural checks on sampled/initial reads
                if checked_pairs < sample_records_to_check:
                    if not (h1.startswith("@") and h2.startswith("@")):
                        raise ValueError(f"Corrupt FASTQ header at record {total_pairs}: R1='{h1.strip()}', R2='{h2.strip()}'")
                    if not (p_plus1.startswith("+") and p_plus2.startswith("+")):
                        raise ValueError(f"Corrupt FASTQ '+' delimiter at record {total_pairs}")
                    if len(s1.strip()) != len(q1.strip()) or len(s2.strip()) != len(q2.strip()):
                        raise ValueError(f"Sequence and quality length mismatch at record {total_pairs}")

                    # Header pairing check: strip /1, /2 or 1:N, 2:N
                    id1 = h1.strip().split()[0].removesuffix("/1")
                    id2 = h2.strip().split()[0].removesuffix("/2")
                    if id1 != id2:
                        raise ValueError(f"Desynchronized read headers at record {total_pairs}: R1='{id1}' != R2='{id2}'")

                    checked_pairs += 1

    except gzip.BadGzipFile as e:
        raise ValueError(f"Corrupted gzip stream in {p1.name} or {p2.name}: {e}")

    return {
        "status": "VALIDATION_PASSED",
        "fastq_1": str(p1),
        "fastq_2": str(p2),
        "total_paired_reads": total_pairs,
        "records_structurally_verified": checked_pairs,
        "is_synchronized": True,
        "read_layout": "paired_end",
        "is_compressed": str(p1).lower().endswith(".gz"),
    }
