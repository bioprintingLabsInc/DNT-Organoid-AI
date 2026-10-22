"""Governed adapter converting native featureCounts output into frozen pipeline matrix."""

import csv
import hashlib
import re
from pathlib import Path
from typing import Any

from .manifest import APPROVED_RUN_ACCESSIONS


def sha256_file(path: str | Path) -> str:
    """Calculate SHA-256 hash of a file."""
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def map_bam_header_to_run_id(col_header: str, approved_runs: tuple[str, ...]) -> str:
    """Deterministically map a BAM column header to one approved SRA Run accession.
    
    Rejects ambiguous or non-matching headers.
    """
    col_str = str(col_header).strip()
    matched = [run for run in approved_runs if re.search(rf"(^|[^A-Za-z0-9]){re.escape(run)}([^A-Za-z0-9]|$)", col_str)]
    if len(matched) == 1:
        return matched[0]
    elif len(matched) > 1:
        raise ValueError(f"Ambiguous BAM column header '{col_header}' matches multiple runs: {matched}")
    else:
        raise ValueError(f"BAM column header '{col_header}' does not match any approved SRA Run accession in {approved_runs}")


def convert_featurecounts_to_matrix(
    raw_featurecounts_path: str | Path,
    output_matrix_path: str | Path,
    output_universe_path: str | Path | None = None,
    approved_runs: tuple[str, ...] = APPROVED_RUN_ACCESSIONS,
) -> dict[str, Any]:
    """Deterministically convert native featureCounts tabular output to frozen-pipeline intake matrix.

    Guarantees:
    - Native file remains untouched as immutable evidence artifact.
    - Strips Chr, Start, End, Strand, Length columns.
    - Maps BAM headers to canonical SRA Run IDs.
    - Preserves exact row count and unique gene identifiers.
    - Asserts 100% numerical conservation (discrepancy count == 0).
    - Verifies pure non-negative integers (no floats, no NaNs).
    """
    raw_path = Path(raw_featurecounts_path)
    out_matrix = Path(output_matrix_path)
    out_universe = Path(output_universe_path) if output_universe_path else None

    if not raw_path.exists():
        raise FileNotFoundError(f"Native featureCounts file not found: {raw_path}")

    raw_sha = sha256_file(raw_path)

    # 1. Read native featureCounts output
    with raw_path.open("r", encoding="utf-8") as f:
        raw_lines = f.readlines()

    if not raw_lines:
        raise ValueError(f"featureCounts output file is empty: {raw_path}")

    # Skip comment lines (typically line 1: # Program:featureCounts ...)
    header_idx = 0
    while header_idx < len(raw_lines) and raw_lines[header_idx].startswith("#"):
        header_idx += 1

    if header_idx >= len(raw_lines):
        raise ValueError(f"No tabular header found in featureCounts file: {raw_path}")

    header_cols = [c.strip() for c in raw_lines[header_idx].strip().split("\t")]

    expected_annot_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    if header_cols[:6] != expected_annot_cols:
        raise ValueError(
            f"Malformed featureCounts header: expected first 6 columns {expected_annot_cols}, got {header_cols[:6]}"
        )

    raw_sample_cols = header_cols[6:]
    if len(raw_sample_cols) != len(approved_runs):
        raise ValueError(
            f"Expected {len(approved_runs)} sample columns, got {len(raw_sample_cols)} in header: {raw_sample_cols}"
        )

    # 2. Map sample columns to approved SRA Run IDs
    mapped_runs = [map_bam_header_to_run_id(col, approved_runs) for col in raw_sample_cols]
    if len(set(mapped_runs)) != len(approved_runs):
        duplicates = [r for r in set(mapped_runs) if mapped_runs.count(r) > 1]
        raise ValueError(f"Duplicate sample mappings detected: {duplicates}")

    # Ensure mapped runs contain all approved runs
    missing_runs = set(approved_runs) - set(mapped_runs)
    if missing_runs:
        raise ValueError(f"Missing approved runs in featureCounts columns: {missing_runs}")

    # 3. Parse data rows, validate values, and perform numerical conservation checks
    gene_ids = []
    matrix_rows = []
    universe_rows = []
    discrepancy_count = 0

    for line_num, line in enumerate(raw_lines[header_idx + 1 :], start=header_idx + 2):
        line_str = line.strip()
        if not line_str:
            continue
        parts = line_str.split("\t")
        if len(parts) != len(header_cols):
            raise ValueError(f"Row width mismatch at line {line_num}: expected {len(header_cols)} cols, got {len(parts)}")

        gid = parts[0].strip()
        if not gid:
            raise ValueError(f"Empty Geneid at line {line_num}")
        gene_ids.append(gid)

        chr_val = parts[1].strip()
        start_val = parts[2].strip()
        end_val = parts[3].strip()
        strand_val = parts[4].strip()
        len_val = parts[5].strip()

        # Parse counts
        sample_counts = []
        row_total_counts = 0
        for s_idx, text_count in enumerate(parts[6:]):
            text_strip = text_count.strip()
            try:
                num = float(text_strip)
            except ValueError:
                raise ValueError(f"Non-numeric count '{text_strip}' for gene '{gid}', column '{mapped_runs[s_idx]}' at line {line_num}")

            if not num.is_integer():
                raise ValueError(f"Non-integer count '{text_strip}' for gene '{gid}', column '{mapped_runs[s_idx]}' at line {line_num}")

            int_val = int(num)
            if int_val < 0:
                raise ValueError(f"Negative count '{int_val}' for gene '{gid}', column '{mapped_runs[s_idx]}' at line {line_num}")

            # Numerical conservation check: string/float to integer exact equality
            if int(text_strip) != int_val:
                discrepancy_count += 1

            sample_counts.append(int_val)
            row_total_counts += int_val

        # Create row dict mapping Run ID -> count
        run_count_dict = dict(zip(mapped_runs, sample_counts))
        # Order counts according to approved_runs order
        ordered_counts = [run_count_dict[run] for run in approved_runs]
        matrix_rows.append((gid, ordered_counts))

        if out_universe:
            universe_rows.append((gid, chr_val, start_val, end_val, strand_val, len_val, row_total_counts))

    # 4. Check uniqueness of Gene IDs
    if len(set(gene_ids)) != len(gene_ids):
        dups = [g for g in set(gene_ids) if gene_ids.count(g) > 1]
        raise ValueError(f"Duplicate gene identifiers found in featureCounts output ({len(dups)} duplicates, e.g. {dups[:5]})")

    # 5. Assert zero discrepancies
    if discrepancy_count != 0:
        raise ValueError(f"Numerical conservation failed with {discrepancy_count} discrepancies between raw and converted counts!")

    # 6. Write frozen-pipeline intake matrix
    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    with out_matrix.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id"] + list(approved_runs))
        for gid, counts in matrix_rows:
            writer.writerow([gid] + counts)

    matrix_sha = sha256_file(out_matrix)

    # 7. Write gene-universe reference table if requested
    universe_sha = ""
    if out_universe:
        out_universe.parent.mkdir(parents=True, exist_ok=True)
        with out_universe.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id", "chr", "start", "end", "strand", "length", "total_counts"])
            for row in universe_rows:
                writer.writerow(row)
        universe_sha = sha256_file(out_universe)

    return {
        "status": "CONVERSION_SUCCESS",
        "raw_featurecounts_path": str(raw_path),
        "raw_featurecounts_sha256": raw_sha,
        "output_matrix_path": str(out_matrix),
        "output_matrix_sha256": matrix_sha,
        "output_universe_path": str(out_universe) if out_universe else None,
        "output_universe_sha256": universe_sha,
        "total_genes_emitted": len(matrix_rows),
        "sample_columns_emitted": list(approved_runs),
        "numerical_conservation_discrepancies": discrepancy_count,
        "is_numerically_conserved": True,
    }
