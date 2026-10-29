"""Deterministic, bounded format detection from file and metadata evidence."""

import csv
import gzip
import io
from pathlib import Path

from src.ingestion.inspection import MAX_INSPECTION_BYTES
from src.ingestion.models import AssetFinding, Confidence, Detection, Inspection, RegisteredAsset


HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
SRA_MAGIC = b"NCBI.sra"

GENE_ALIASES = {"gene", "gene_id", "gene_name", "gene_symbol", "symbol", "ensembl_id", "id", "feature_id"}
EFFECT_ALIASES = {"log2foldchange", "logfc", "log2fc", "lfc", "fold_change", "foldchange", "log2_fc", "beta"}
STAT_ALIASES = {"pvalue", "p_value", "p.value", "pval", "p_val", "padj", "p_adj", "adj.p.val", "fdr", "qval", "qvalue", "q_value"}


def detect_format(asset: RegisteredAsset, inspection: Inspection) -> Detection:
    """Detect supported formats conservatively without loading biological data."""
    if inspection.findings:
        return Detection(None, Confidence.NONE, (), inspection.findings)
    suffixes = inspection.suffixes
    if _is_fastq_name(suffixes):
        return _detect_fastq(asset.local_path, inspection)
    if suffixes[-1:] == (".sra",):
        if inspection.prefix.startswith(SRA_MAGIC):
            return Detection("SRA", Confidence.HIGH, (".sra extension", "NCBI SRA signature"))
        return _review("SRA extension lacks a recognized local SRA signature", "sra_unconfirmed")
    if suffixes[-1:] == (".h5ad",):
        if inspection.prefix.startswith(HDF5_MAGIC):
            return Detection("H5AD", Confidence.HIGH, (".h5ad extension", "HDF5 signature"))
        return _review(".h5ad extension lacks an HDF5 signature", "h5ad_unconfirmed")
    if suffixes[-1:] == (".rds",):
        if inspection.prefix.startswith((b"X\n", b"A\n", b"B\n")) or inspection.prefix.startswith(b"\x1f\x8b"):
            return Detection("RDS", Confidence.MEDIUM, (".rds extension", "R serialization-compatible signature"))
        return _review("RDS contents cannot be identified safely", "rds_unconfirmed")
    if suffixes[-1:] == (".xlsx",):
        return _detect_xlsx(asset, inspection)
    if _is_tabular_name(suffixes):
        return _detect_matrix(asset, inspection)
    return _review("No supported format can be established from bounded evidence", "unsupported_or_unknown")


def _detect_xlsx(asset: RegisteredAsset, inspection: Inspection) -> Detection:
    from src.ingestion.xlsx_reader import classify_xlsx_workbook
    selected, all_sheets, findings, prov = classify_xlsx_workbook(asset.local_path)
    if findings:
        return Detection(None, Confidence.LOW, (), tuple(findings))

    evidence = (
        ".xlsx workbook",
        f"selected count sheet '{selected}'",
        f"{len(all_sheets)} total sheets",
        "non-negative integer-like values",
    )
    if asset.declared_format == "raw_counts" and asset.modalities == ("bulk_rna_seq",):
        return Detection("raw_counts", Confidence.HIGH, evidence + ("declared raw_counts with bulk modality",))
    return Detection("raw_counts", Confidence.MEDIUM, evidence)


def _detect_fastq(path: Path, inspection: Inspection) -> Detection:
    try:
        if inspection.suffixes[-1:] == (".gz",):
            with gzip.open(path, "rb") as stream:
                content = stream.read(MAX_INSPECTION_BYTES + 1)
            if len(content) > MAX_INSPECTION_BYTES:
                content = content[:MAX_INSPECTION_BYTES]
        else:
            content = inspection.prefix
        lines = content.decode("utf-8", errors="strict").splitlines()[:8]
    except (OSError, UnicodeError, EOFError) as error:
        return Detection("FASTQ", Confidence.LOW, ("FASTQ filename extension",), (
            AssetFinding("ERROR", "malformed_fastq", "source_reference", f"FASTQ cannot be decoded: {error}"),
        ))
    records = [lines[index:index + 4] for index in range(0, len(lines), 4) if lines[index]]
    valid = bool(records) and all(len(record) == 4 and record[0].startswith("@") and record[2].startswith("+") and bool(record[1]) and len(record[1]) == len(record[3]) for record in records)
    if not valid:
        return Detection("FASTQ", Confidence.LOW, ("FASTQ filename extension",), (
            AssetFinding("ERROR", "malformed_fastq", "source_reference", "FASTQ structure is malformed in the inspected records."),
        ))
    compression = "gzip compression" if inspection.suffixes[-1:] == (".gz",) else "uncompressed text"
    return Detection("FASTQ", Confidence.HIGH, ("FASTQ filename extension", compression, "valid four-line FASTQ records"))


def _detect_matrix(asset: RegisteredAsset, inspection: Inspection) -> Detection:
    try:
        is_gz = inspection.suffixes[-1:] == (".gz",)
        if is_gz:
            with gzip.open(asset.local_path, "rb") as stream:
                content = stream.read(MAX_INSPECTION_BYTES + 1)
            if len(content) > MAX_INSPECTION_BYTES:
                content = content[:MAX_INSPECTION_BYTES]
            text = content.decode("utf-8-sig", errors="replace")
            delimiter = "," if len(inspection.suffixes) >= 2 and inspection.suffixes[-2] == ".csv" else "\t"
        else:
            text = inspection.prefix.decode("utf-8-sig")
            delimiter = "\t" if inspection.suffixes[-1] in {".tsv", ".txt"} else ","
        rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))[:51]
    except (UnicodeError, csv.Error, OSError) as error:
        return _review(f"Tabular file cannot be inspected: {error}", "malformed_matrix")
    if len(rows) < 2 or len(rows[0]) < 2:
        return _review("Matrix requires a header and at least one data row", "malformed_matrix")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        return _review("Matrix rows have inconsistent column counts", "malformed_matrix")

    # Check for differential expression summary table structure
    deg_result = _check_deg_table(rows)
    if deg_result is not None:
        return deg_result

    try:
        values = [float(value) for row in rows[1:] for value in row[1:]]
    except ValueError:
        return _review("Matrix body is not consistently numeric", "ambiguous_matrix")
    if not values:
        return _review("Matrix has no numeric expression values", "empty_matrix")
    integer_like = all(value >= 0 and value.is_integer() for value in values)
    compression_ev = "gzip-compressed tabular structure" if is_gz else "tabular matrix structure"
    evidence = (compression_ev, f"{len(rows) - 1} inspected data rows", "non-negative integer-like values" if integer_like else "non-integer numeric values")
    if not integer_like:
        return Detection("processed_expression_matrix", Confidence.HIGH, evidence)
    if asset.declared_format == "raw_counts" and asset.modalities == ("bulk_rna_seq",):
        return Detection("raw_counts", Confidence.MEDIUM, evidence + ("declared raw_counts with bulk modality",))
    if asset.declared_format == "UMI_counts" and asset.modalities and set(asset.modalities) <= {"scRNA_seq", "snRNA_seq"}:
        return Detection("UMI_counts", Confidence.MEDIUM, evidence + ("declared UMI_counts with single-cell/nucleus modality",))
    return Detection(None, Confidence.LOW, evidence, (
        AssetFinding("WARNING", "ambiguous_integer_matrix", "source_reference", "Integer-like values do not distinguish raw_counts from UMI_counts without compatible metadata."),
    ))


def _check_deg_table(rows: list[list[str]]) -> Detection | None:
    """Classify differential expression summary tables by recognizing key column signatures."""
    if len(rows) < 2 or len(rows[0]) < 3:
        return None
    header = [h.strip().lower() for h in rows[0]]
    gene_idx = next((i for i, h in enumerate(header) if h in GENE_ALIASES), None)
    effect_idx = next((i for i, h in enumerate(header) if h in EFFECT_ALIASES), None)
    stat_idx = next((i for i, h in enumerate(header) if h in STAT_ALIASES), None)

    if gene_idx is None or effect_idx is None or stat_idx is None:
        return None

    valid_rows = 0
    for row in rows[1:]:
        if len(row) <= max(gene_idx, effect_idx, stat_idx):
            continue
        try:
            float(row[effect_idx].strip())
            float(row[stat_idx].strip())
            valid_rows += 1
        except ValueError:
            continue

    if valid_rows == 0:
        return None

    evidence = (
        "tabular differential expression table structure",
        f"{len(rows) - 1} inspected rows",
        f"identified gene identifier column '{rows[0][gene_idx]}'",
        f"identified effect size column '{rows[0][effect_idx]}'",
        f"identified significance column '{rows[0][stat_idx]}'",
    )
    return Detection("differential_expression_table", Confidence.HIGH, evidence)


def _is_fastq_name(suffixes: tuple[str, ...]) -> bool:
    return suffixes[-1:] in ((".fastq",), (".fq",)) or suffixes[-2:] in ((".fastq", ".gz"), (".fq", ".gz"))


def _is_tabular_name(suffixes: tuple[str, ...]) -> bool:
    if not suffixes:
        return False
    if suffixes[-1] in {".csv", ".tsv", ".txt"}:
        return True
    if len(suffixes) >= 2 and suffixes[-1] == ".gz" and suffixes[-2] in {".csv", ".tsv", ".txt"}:
        return True
    return False


def _review(message: str, rule: str) -> Detection:
    return Detection(None, Confidence.LOW, (), (AssetFinding("WARNING", rule, "source_reference", message),))
