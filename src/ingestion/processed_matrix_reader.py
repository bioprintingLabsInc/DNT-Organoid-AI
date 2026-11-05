"""Safe, bounded, and deterministic ingestion of processed expression matrices.

Implements Phase U3.2 processed expression matrix intake:
1. Supported tabular formats (.tsv, .csv, .txt, .tsv.gz, .csv.gz, and .xlsx).
2. Exact preservation of floating-point numerical values, gene identifiers, and sample headers.
3. Zero data fabrication: no integer conversion, rounding, or silent imputation.
4. Structural validation: matrix dimensions, missing values, non-numeric values, and duplicate identifiers.
5. Deterministic row ordering and column representation.
"""

from collections import Counter
import csv
from dataclasses import dataclass, field
import gzip
import io
import math
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

from src.ingestion.models import AssetFinding


@dataclass(frozen=True, slots=True)
class ProcessedExpressionMatrix:
    """Continuous sample-level processed bulk RNA-seq expression matrix.

    Preserves exact floating-point numerical values, original gene IDs,
    and original sample identifiers without rounding or integer conversion.
    """

    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    columns: tuple[tuple[float, ...], ...]  # columns[sample_index][gene_index]

    @property
    def number_of_genes(self) -> int:
        return len(self.gene_ids)

    @property
    def number_of_samples(self) -> int:
        return len(self.sample_ids)

    def get_value(self, gene_id: str, sample_id: str) -> float | None:
        """Retrieve continuous expression value for a given gene and sample."""
        try:
            s_idx = self.sample_ids.index(sample_id)
            g_idx = self.gene_ids.index(gene_id)
            return self.columns[s_idx][g_idx]
        except ValueError:
            return None

    def get_sample_column(self, sample_id: str) -> tuple[float, ...] | None:
        """Retrieve all gene expression values for a specified sample."""
        try:
            s_idx = self.sample_ids.index(sample_id)
            return self.columns[s_idx]
        except ValueError:
            return None

    def get_gene_row(self, gene_id: str) -> tuple[float, ...] | None:
        """Retrieve expression values across all samples for a specified gene."""
        try:
            g_idx = self.gene_ids.index(gene_id)
            return tuple(col[g_idx] for col in self.columns)
        except ValueError:
            return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "number_of_genes": self.number_of_genes,
            "number_of_samples": self.number_of_samples,
            "gene_ids": list(self.gene_ids),
            "sample_ids": list(self.sample_ids),
        }


@dataclass(frozen=True, slots=True)
class ProcessedMatrixData:
    """Container for parsed processed expression matrix and intake audit findings."""

    matrix: ProcessedExpressionMatrix | None
    findings: tuple[AssetFinding, ...]
    source_path: Path
    sheet_name: str | None = None
    raw_row_count: int = 0
    raw_column_count: int = 0

    @property
    def is_valid(self) -> bool:
        """Valid if matrix exists and no ERROR-level findings were recorded."""
        return self.matrix is not None and not any(f.severity == "ERROR" for f in self.findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_valid": self.is_valid,
            "source_path": str(self.source_path),
            "sheet_name": self.sheet_name,
            "raw_row_count": self.raw_row_count,
            "raw_column_count": self.raw_column_count,
            "findings": [
                {
                    "severity": f.severity,
                    "rule": f.rule,
                    "path": f.path,
                    "message": f.message,
                }
                for f in self.findings
            ],
            "matrix": self.matrix.to_dict() if self.matrix else None,
        }


def read_processed_matrix(
    path: str | Path,
    sheet_name: str | None = None,
) -> ProcessedMatrixData:
    """Read and validate a processed bulk RNA-seq expression matrix from disk.

    Supports .tsv, .csv, .txt, .tsv.gz, .csv.gz, and .xlsx files.
    Preserves exact floating-point numerical values, gene identifiers, and sample headers.
    """
    path_p = Path(path)
    if not path_p.is_file():
        finding = AssetFinding("ERROR", "file_not_found", "source_reference", f"File not found: {path_p}")
        return ProcessedMatrixData(None, (finding,), path_p, sheet_name)

    suffix_lower = path_p.name.lower()
    if suffix_lower.endswith(".xlsx"):
        return _read_processed_xlsx(path_p, sheet_name=sheet_name)
    else:
        return _read_processed_tabular(path_p)


def _read_processed_tabular(path: Path) -> ProcessedMatrixData:
    """Read delimited text matrix (CSV/TSV/TXT/gz)."""
    findings: list[AssetFinding] = []
    is_gz = path.name.lower().endswith(".gz")

    try:
        if is_gz:
            with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as f:
                content = f.read()
        else:
            with open(path, "rt", encoding="utf-8-sig", errors="replace") as f:
                content = f.read()
    except Exception as e:
        finding = AssetFinding("ERROR", "cannot_read_file", "source_reference", f"Failed to read file: {e}")
        return ProcessedMatrixData(None, (finding,), path)

    if not content or not content.strip():
        finding = AssetFinding("ERROR", "empty_file", "source_reference", "File is empty.")
        return ProcessedMatrixData(None, (finding,), path)

    # Determine delimiter
    first_line = content.splitlines()[0] if content.splitlines() else ""
    if path.name.lower().endswith(".csv") or path.name.lower().endswith(".csv.gz"):
        delimiter = ","
    elif "\t" in first_line:
        delimiter = "\t"
    elif "," in first_line:
        delimiter = ","
    else:
        delimiter = "\t"

    try:
        reader = csv.reader(io.StringIO(content), delimiter=delimiter)
        raw_rows = list(reader)
    except Exception as e:
        finding = AssetFinding("ERROR", "csv_parse_error", "source_reference", f"CSV parsing error: {e}")
        return ProcessedMatrixData(None, (finding,), path)

    return _parse_rows_into_matrix(raw_rows, path)


def _read_processed_xlsx(path: Path, sheet_name: str | None = None) -> ProcessedMatrixData:
    """Read Excel workbook containing a processed expression matrix."""
    from src.ingestion.xlsx_reader import NS_MAIN, inspect_xlsx_workbook, parse_cell_ref, read_shared_strings

    findings: list[AssetFinding] = []
    try:
        is_xlsx, _, sheets = inspect_xlsx_workbook(path)
        if not is_xlsx or not sheets:
            finding = AssetFinding("ERROR", "malformed_xlsx", "source_reference", "File is not a valid XLSX workbook.")
            return ProcessedMatrixData(None, (finding,), path)

        target_sheet = None
        if sheet_name:
            target_sheet = next((s for s in sheets if s["name"].strip().lower() == sheet_name.strip().lower()), None)
            if not target_sheet:
                finding = AssetFinding(
                    "ERROR",
                    "sheet_not_found",
                    "source_reference",
                    f"Requested sheet '{sheet_name}' not found in workbook. Available: {[s['name'] for s in sheets]}",
                )
                return ProcessedMatrixData(None, (finding,), path)
        else:
            target_sheet = sheets[0]

        with zipfile.ZipFile(path) as z:
            shared_strings = read_shared_strings(z)
            sheet_xml = z.read(target_sheet["target"])
            root = ET.fromstring(sheet_xml)

            rows: list[list[str]] = []
            row_elements = root.findall(f".//{{{NS_MAIN}}}row")
            for r_elem in row_elements:
                cell_elements = r_elem.findall(f"{{{NS_MAIN}}}c")
                row_dict: dict[int, str] = {}
                for c_elem in cell_elements:
                    ref = c_elem.get("r", "")
                    col_idx, _ = parse_cell_ref(ref)
                    t_attr = c_elem.get("t", "")
                    v_elem = c_elem.find(f"{{{NS_MAIN}}}v")
                    val = v_elem.text if v_elem is not None else ""
                    if t_attr == "s" and val.isdigit() and int(val) < len(shared_strings):
                        val = shared_strings[int(val)]
                    row_dict[col_idx] = val.strip()

                if row_dict:
                    max_col = max(row_dict.keys())
                    dense_row = [row_dict.get(c, "") for c in range(1, max_col + 1)]
                    rows.append(dense_row)

        return _parse_rows_into_matrix(rows, path, sheet_name=target_sheet["name"])

    except Exception as e:
        finding = AssetFinding("ERROR", "xlsx_read_error", "source_reference", f"Failed to read Excel workbook: {e}")
        return ProcessedMatrixData(None, (finding,), path, sheet_name)


def _parse_rows_into_matrix(
    rows: list[list[str]],
    path: Path,
    sheet_name: str | None = None,
) -> ProcessedMatrixData:
    """Parse raw row lists into a validated ProcessedExpressionMatrix with full audit findings."""
    findings: list[AssetFinding] = []
    if not rows or not any(rows):
        finding = AssetFinding("ERROR", "empty_matrix", "matrix", "Matrix contains no data rows.")
        return ProcessedMatrixData(None, (finding,), path, sheet_name)

    # Filter leading empty rows
    start_idx = 0
    while start_idx < len(rows) and not any(cell.strip() for cell in rows[start_idx]):
        start_idx += 1

    if start_idx >= len(rows):
        finding = AssetFinding("ERROR", "empty_matrix", "matrix", "Matrix contains only whitespace rows.")
        return ProcessedMatrixData(None, (finding,), path, sheet_name)

    header = [cell.strip() for cell in rows[start_idx]]
    if not header or not header[0]:
        findings.append(
            AssetFinding("ERROR", "gene_identifier_column_required", "header[0]", "Gene identifier column is missing.")
        )

    samples = tuple(header[1:])
    if not samples:
        findings.append(
            AssetFinding("ERROR", "sample_columns_required", "header", "At least one sample column is required.")
        )

    sample_counts = Counter(samples)
    for s in sorted([x for x, c in sample_counts.items() if c > 1]):
        findings.append(
            AssetFinding("ERROR", "unique_sample_identifiers_required", f"header.{s}", f"Duplicate sample identifier '{s}'.")
        )

    gene_ids: list[str] = []
    values_by_sample: list[list[float]] = [[] for _ in samples]
    raw_data_rows = 0

    for row_idx, row in enumerate(rows[start_idx + 1 :], start=start_idx + 2):
        if not row or not any(cell.strip() for cell in row):
            findings.append(
                AssetFinding("ERROR", "empty_row", f"row[{row_idx}]", f"Empty row at index {row_idx}.")
            )
            continue

        raw_data_rows += 1

        if len(row) != len(header):
            findings.append(
                AssetFinding(
                    "ERROR",
                    "malformed_row",
                    f"row[{row_idx}]",
                    f"Row {row_idx} width ({len(row)}) differs from header width ({len(header)}).",
                )
            )
            continue

        gene_id = row[0].strip()
        if not gene_id:
            findings.append(
                AssetFinding("ERROR", "empty_gene_identifier", f"row[{row_idx}].gene_id", f"Empty gene identifier at row {row_idx}.")
            )

        gene_ids.append(gene_id)

        for s_idx, text in enumerate(row[1:]):
            sample_name = samples[s_idx]
            cell_str = text.strip()
            cell_path = f"row[{row_idx}].{sample_name}"

            if cell_str == "" or cell_str.lower() in {"na", "nan", "null", "none"}:
                findings.append(
                    AssetFinding("ERROR", "missing_expression_values_forbidden", cell_path, f"Missing expression value for gene '{gene_id}' in sample '{sample_name}'.")
                )
                continue

            try:
                num = float(cell_str)
            except ValueError:
                findings.append(
                    AssetFinding("ERROR", "numeric_expression_required", cell_path, f"Expression value '{cell_str}' is not numeric.")
                )
                continue

            if math.isnan(num) or math.isinf(num):
                findings.append(
                    AssetFinding("ERROR", "non_finite_expression_values_forbidden", cell_path, f"Expression value '{cell_str}' is non-finite.")
                )
                continue

            values_by_sample[s_idx].append(num)

    # Check for duplicate gene identifiers (REVIEW severity - preserved without silent deduplication)
    gene_counts = Counter(gene_ids)
    for g in sorted([x for x, c in gene_counts.items() if c > 1]):
        findings.append(
            AssetFinding(
                "REVIEW",
                "duplicate_gene_identifiers",
                f"gene[{g}]",
                f"Duplicate gene identifier '{g}' detected ({gene_counts[g]} occurrences); all rows preserved without aggregation.",
            )
        )

    # If any ERROR was recorded, matrix cannot be constructed
    if any(f.severity == "ERROR" for f in findings):
        return ProcessedMatrixData(
            None,
            tuple(findings),
            path,
            sheet_name,
            raw_row_count=raw_data_rows,
            raw_column_count=len(header),
        )

    matrix = ProcessedExpressionMatrix(
        gene_ids=tuple(gene_ids),
        sample_ids=samples,
        columns=tuple(tuple(col) for col in values_by_sample),
    )

    return ProcessedMatrixData(
        matrix=matrix,
        findings=tuple(findings),
        source_path=path,
        sheet_name=sheet_name,
        raw_row_count=raw_data_rows,
        raw_column_count=len(header),
    )
