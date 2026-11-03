"""Safe, bounded, and deterministic ingestion of published differential expression summary tables.

Implements Phase U3.1 ingestion:
1. Tabular (CSV, TSV, TXT, optionally gzip-compressed) and Excel (.xlsx) parsing.
2. Conservative column header alias discovery without guessing.
3. Exact preservation of published numeric values (log2FC, p-value, padj/FDR, SE, stat, baseMean).
4. Strictly zero data fabrication: missing standard errors, test statistics, base means,
   or omitted/NA p-values remain None.
5. Preserves full provenance, source row ordering, and row index mapping.
"""

import csv
from dataclasses import dataclass, field
import gzip
import io
import math
from pathlib import Path
import re
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

from src.ingestion.models import AssetFinding
from src.ingestion.xlsx_reader import NS_MAIN, inspect_xlsx_workbook, parse_cell_ref, read_shared_strings

# Alias vocabularies for column discovery
GENE_ID_ALIASES = {
    "gene_id", "geneid", "gene", "ensembl_id", "ensembl_gene_id", "ensembl",
    "feature_id", "id", "transcript_id", "probeset_id", "probe_id", "target_id",
}
GENE_SYMBOL_ALIASES = {
    "gene_name", "genename", "gene_symbol", "symbol", "hgnc_symbol", "approved_symbol",
}
EFFECT_ALIASES = {
    "log2foldchange", "log2_fold_change", "logfc", "log2fc", "lfc",
    "fold_change", "foldchange", "log2_fc", "beta", "estimate", "log2(fold_change)",
}
PVAL_ALIASES = {
    "pvalue", "p_value", "p.value", "pval", "p_val", "p",
}
PADJ_ALIASES = {
    "padj", "p_adj", "adj.p.val", "adj_p_val", "adjpval", "fdr",
    "qval", "qvalue", "q_value", "p_adjusted", "adj_p_value", "p.adjust",
}
SE_ALIASES = {
    "lfcse", "lfc_se", "se", "stderr", "standard_error", "sem",
}
STAT_ALIASES = {
    "stat", "wald_statistic", "wald_stat", "t", "z", "f", "lr", "chisq",
}
BASE_MEAN_ALIASES = {
    "basemean", "base_mean", "aveexpr", "ave_expr", "mean_expression", "avg_expr", "mean_expr", "logcpm",
}


@dataclass(frozen=True, slots=True)
class DegTableRow:
    """A single row from a published differential expression table."""

    source_index: int  # 1-based source row index in the table
    raw_gene_id: str
    raw_gene_symbol: str | None
    log2_fold_change: float | None
    p_value: float | None
    adjusted_p_value: float | None
    lfc_standard_error: float | None = None
    wald_statistic: float | None = None
    base_mean: float | None = None
    extra_fields: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class DegTableData:
    """Parsed and validated differential expression table."""

    rows: tuple[DegTableRow, ...]
    total_rows: int
    column_mapping: dict[str, str]  # role -> source column header name
    source_path: Path
    sheet_name: str | None
    findings: tuple[AssetFinding, ...]
    provenance: dict[str, Any]

    @property
    def is_valid(self) -> bool:
        """True if no fatal errors were found and at least one data row was extracted."""
        return not any(f.severity == "ERROR" for f in self.findings) and len(self.rows) > 0


def _clean_header_key(header: str) -> str:
    """Normalize a header name for robust alias comparison."""
    return re.sub(r"[^a-z0-9_.]", "", header.strip().lower())


def discover_deg_columns(
    headers: list[str],
    overrides: dict[str, str] | None = None,
) -> tuple[dict[str, int], dict[str, str], list[AssetFinding]]:
    """Identify column roles from header names or explicit caller overrides.

    Returns:
        (col_indices, role_to_header_name, findings)
    """
    findings: list[AssetFinding] = []
    col_indices: dict[str, int] = {}
    col_mapping: dict[str, str] = {}

    cleaned_headers = [_clean_header_key(h) for h in headers]
    hdr_map = {clean: idx for idx, clean in enumerate(cleaned_headers)}

    # Apply caller overrides if provided
    if overrides:
        for role, name in overrides.items():
            c_name = _clean_header_key(name)
            if c_name in hdr_map:
                col_indices[role] = hdr_map[c_name]
                col_mapping[role] = headers[hdr_map[c_name]]
            else:
                findings.append(
                    AssetFinding(
                        "ERROR",
                        "override_column_not_found",
                        "deg_table",
                        f"Explicit override column '{name}' for role '{role}' not found in headers.",
                    )
                )

    # 1. Gene Identifier
    if "gene_id" not in col_indices:
        # Check specific ID aliases first
        for idx, h in enumerate(cleaned_headers):
            if h in GENE_ID_ALIASES:
                col_indices["gene_id"] = idx
                col_mapping["gene_id"] = headers[idx]
                break

    # 2. Gene Symbol (separate from gene_id if both exist)
    if "gene_symbol" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in GENE_SYMBOL_ALIASES and idx != col_indices.get("gene_id"):
                col_indices["gene_symbol"] = idx
                col_mapping["gene_symbol"] = headers[idx]
                break

    # If no gene_id found yet, but a gene_symbol column exists, use it as primary identifier
    if "gene_id" not in col_indices and "gene_symbol" in col_indices:
        col_indices["gene_id"] = col_indices["gene_symbol"]
        col_mapping["gene_id"] = col_mapping["gene_symbol"]

    # If still no gene identifier column, check if any symbol alias matches
    if "gene_id" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in GENE_SYMBOL_ALIASES:
                col_indices["gene_id"] = idx
                col_mapping["gene_id"] = headers[idx]
                break

    if "gene_id" not in col_indices:
        findings.append(
            AssetFinding(
                "ERROR",
                "missing_gene_column",
                "deg_table",
                "No recognizable gene identifier column found (e.g. gene_id, gene_name, symbol, ensembl_id).",
            )
        )

    # 3. Effect Size (log2 fold change)
    if "log2fc" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in EFFECT_ALIASES:
                col_indices["log2fc"] = idx
                col_mapping["log2fc"] = headers[idx]
                break

    if "log2fc" not in col_indices:
        findings.append(
            AssetFinding(
                "ERROR",
                "missing_effect_size_column",
                "deg_table",
                "No recognizable log2 fold change or effect size column found (e.g. log2FoldChange, logFC, lfc).",
            )
        )

    # 4. Significance: Adjusted P-value / FDR
    if "padj" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in PADJ_ALIASES:
                col_indices["padj"] = idx
                col_mapping["padj"] = headers[idx]
                break

    # 5. Significance: Nominal P-value
    if "pvalue" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in PVAL_ALIASES and idx != col_indices.get("padj"):
                col_indices["pvalue"] = idx
                col_mapping["pvalue"] = headers[idx]
                break

    # If neither pvalue nor padj found, issue a warning (some papers publish unthresholded effect size tables)
    if "padj" not in col_indices and "pvalue" not in col_indices:
        findings.append(
            AssetFinding(
                "WARNING",
                "missing_significance_columns",
                "deg_table",
                "Neither nominal p-value nor adjusted p-value/FDR column was discovered.",
            )
        )

    # 6. Standard Error (Optional)
    if "lfc_se" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in SE_ALIASES:
                col_indices["lfc_se"] = idx
                col_mapping["lfc_se"] = headers[idx]
                break

    # 7. Test Statistic (Optional)
    if "stat" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in STAT_ALIASES:
                col_indices["stat"] = idx
                col_mapping["stat"] = headers[idx]
                break

    # 8. Base Mean / Average Expression (Optional)
    if "base_mean" not in col_indices:
        for idx, h in enumerate(cleaned_headers):
            if h in BASE_MEAN_ALIASES:
                col_indices["base_mean"] = idx
                col_mapping["base_mean"] = headers[idx]
                break

    return col_indices, col_mapping, findings


def _parse_numeric(val: Any) -> tuple[float | None, bool]:
    """Parse numeric float from cell string or primitive.

    Returns: (parsed_float, is_unparsable_string)
    """
    if val is None:
        return None, False
    if isinstance(val, (int, float)):
        if math.isnan(val):
            return None, False
        return float(val), False

    s = str(val).strip()
    if s == "" or s.lower() in {"na", "nan", "null", "none", ".", "-"}:
        return None, False
    if s.lower() in {"inf", "+inf"}:
        return float("inf"), False
    if s.lower() == "-inf":
        return float("-inf"), False

    try:
        return float(s), False
    except ValueError:
        return None, True


def read_deg_table(
    path: str | Path,
    sheet_name: str | None = None,
    column_mapping: dict[str, str] | None = None,
) -> DegTableData:
    """Read and validate a published differential expression table from disk.

    Supports .tsv, .csv, .txt, .tsv.gz, .csv.gz, and .xlsx files.
    """
    path_p = Path(path)
    if not path_p.is_file():
        finding = AssetFinding("ERROR", "file_not_found", "source_reference", f"File not found: {path_p}")
        return DegTableData((), 0, {}, path_p, sheet_name, (finding,), {"status": "FILE_NOT_FOUND"})

    suffix_lower = path_p.name.lower()
    if suffix_lower.endswith(".xlsx"):
        return _read_deg_xlsx(path_p, sheet_name=sheet_name, overrides=column_mapping)
    else:
        return _read_deg_tabular(path_p, overrides=column_mapping)


def _read_deg_tabular(
    path: Path,
    overrides: dict[str, str] | None = None,
) -> DegTableData:
    """Read a delimited text differential expression table (CSV/TSV/TXT/gz)."""
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
        return DegTableData((), 0, {}, path, None, (finding,), {"status": "READ_ERROR"})

    # Determine delimiter
    lines = [line for line in content.splitlines() if line.strip()]
    if not lines:
        finding = AssetFinding("ERROR", "empty_file", "source_reference", "File contains no non-empty lines.")
        return DegTableData((), 0, {}, path, None, (finding,), {"status": "EMPTY_FILE"})

    first_line = lines[0]
    delimiter = "\t" if "\t" in first_line else ","

    reader = csv.reader(io.StringIO(content), delimiter=delimiter)
    all_rows = list(reader)
    if len(all_rows) < 2:
        finding = AssetFinding("ERROR", "insufficient_rows", "deg_table", "Table must contain a header and at least one data row.")
        return DegTableData((), 0, {}, path, None, (finding,), {"status": "INSUFFICIENT_ROWS"})

    headers = [h.strip() for h in all_rows[0]]
    col_indices, col_mapping, col_findings = discover_deg_columns(headers, overrides=overrides)
    findings.extend(col_findings)

    if any(f.severity == "ERROR" for f in findings):
        return DegTableData((), len(all_rows) - 1, col_mapping, path, None, tuple(findings), {"status": "COLUMN_DISCOVERY_FAILED"})

    parsed_rows, row_findings = _parse_rows_from_list(
        all_rows[1:],
        headers,
        col_indices,
    )
    findings.extend(row_findings)

    provenance = {
        "source_format": "tabular_delimited",
        "is_gzip": is_gz,
        "delimiter": "tab" if delimiter == "\t" else "comma",
        "total_input_rows": len(parsed_rows),
        "column_mapping": col_mapping,
        "path": str(path),
    }

    return DegTableData(
        rows=tuple(parsed_rows),
        total_rows=len(parsed_rows),
        column_mapping=col_mapping,
        source_path=path,
        sheet_name=None,
        findings=tuple(findings),
        provenance=provenance,
    )


def _read_deg_xlsx(
    path: Path,
    sheet_name: str | None = None,
    overrides: dict[str, str] | None = None,
) -> DegTableData:
    """Read a differential expression sheet from an Excel (.xlsx) workbook."""
    findings: list[AssetFinding] = []
    is_valid, err, sheets = inspect_xlsx_workbook(path)
    if not is_valid:
        finding = AssetFinding("ERROR", "corrupt_workbook", "source_reference", f"Invalid XLSX workbook: {err}")
        return DegTableData((), 0, {}, path, sheet_name, (finding,), {"status": "CORRUPT_WORKBOOK"})

    with zipfile.ZipFile(path, "r") as z:
        shared = read_shared_strings(z)

        # 1. Identify target sheet
        target_sheet_info = None
        if sheet_name:
            for s in sheets:
                if s["name"].lower() == sheet_name.strip().lower():
                    target_sheet_info = s
                    break
            if target_sheet_info is None:
                finding = AssetFinding(
                    "ERROR",
                    "sheet_not_found",
                    "source_reference",
                    f"Requested sheet '{sheet_name}' not found in workbook (available: {[s['name'] for s in sheets]}).",
                )
                return DegTableData((), 0, {}, path, sheet_name, (finding,), {"status": "SHEET_NOT_FOUND"})
        else:
            # Enumerate candidate sheets that have DEG headers
            candidate_sheets = []
            for s in sheets:
                header = _peek_xlsx_sheet_header(z, s["target"], shared)
                if header:
                    cols, _, errs = discover_deg_columns(header, overrides=overrides)
                    if not any(e.severity == "ERROR" for e in errs):
                        candidate_sheets.append(s)

            if len(candidate_sheets) == 1:
                target_sheet_info = candidate_sheets[0]
            elif len(candidate_sheets) > 1:
                finding = AssetFinding(
                    "ERROR",
                    "ambiguous_deg_sheets",
                    "source_reference",
                    f"Multiple candidate differential expression sheets found: {[s['name'] for s in candidate_sheets]}. Please specify sheet_name explicitly.",
                )
                return DegTableData((), 0, {}, path, None, (finding,), {"status": "AMBIGUOUS_DEG_SHEETS"})
            else:
                finding = AssetFinding(
                    "ERROR",
                    "no_deg_sheet_found",
                    "source_reference",
                    f"No sheet containing differential expression columns found in workbook (inspected {len(sheets)} sheets).",
                )
                return DegTableData((), 0, {}, path, None, (finding,), {"status": "NO_DEG_SHEET_FOUND"})

        actual_sheet_name = target_sheet_info["name"]

        # 2. Parse selected sheet
        headers, data_rows = _parse_full_xlsx_sheet(z, target_sheet_info["target"], shared)
        if not headers or not data_rows:
            finding = AssetFinding("ERROR", "empty_sheet", "deg_table", f"Sheet '{actual_sheet_name}' contains insufficient data.")
            return DegTableData((), 0, {}, path, actual_sheet_name, (finding,), {"status": "EMPTY_SHEET"})

        col_indices, col_mapping, col_findings = discover_deg_columns(headers, overrides=overrides)
        findings.extend(col_findings)

        if any(f.severity == "ERROR" for f in findings):
            return DegTableData((), len(data_rows), col_mapping, path, actual_sheet_name, tuple(findings), {"status": "COLUMN_DISCOVERY_FAILED"})

        parsed_rows, row_findings = _parse_rows_from_list(
            data_rows,
            headers,
            col_indices,
        )
        findings.extend(row_findings)

        provenance = {
            "source_format": "xlsx_workbook",
            "selected_sheet": actual_sheet_name,
            "all_sheets": [s["name"] for s in sheets],
            "total_input_rows": len(parsed_rows),
            "column_mapping": col_mapping,
            "path": str(path),
        }

        return DegTableData(
            rows=tuple(parsed_rows),
            total_rows=len(parsed_rows),
            column_mapping=col_mapping,
            source_path=path,
            sheet_name=actual_sheet_name,
            findings=tuple(findings),
            provenance=provenance,
        )


def _peek_xlsx_sheet_header(
    z: zipfile.ZipFile,
    target: str,
    shared_strings: list[str],
) -> list[str]:
    """Read only the first row (headers) of a worksheet."""
    if target not in z.namelist():
        return []
    try:
        with z.open(target) as stream:
            for event, elem in ET.iterparse(stream, events=("end",)):
                if elem.tag == f"{{{NS_MAIN}}}row":
                    r_cells = {}
                    for c in elem.findall(f"{{{NS_MAIN}}}c"):
                        ref = c.attrib.get("r", "")
                        c_type = c.attrib.get("t", "")
                        _, col_idx = parse_cell_ref(ref)
                        v_elem = c.find(f"{{{NS_MAIN}}}v")
                        raw_val = v_elem.text if v_elem is not None else None
                        cell_val = ""
                        if c_type == "s" and raw_val is not None:
                            s_idx = int(raw_val)
                            cell_val = shared_strings[s_idx] if s_idx < len(shared_strings) else ""
                        elif c_type == "inlineStr":
                            is_el = c.find(f"{{{NS_MAIN}}}is")
                            if is_el is not None:
                                t_el = is_el.find(f"{{{NS_MAIN}}}t")
                                if t_el is not None and t_el.text:
                                    cell_val = t_el.text
                        elif raw_val is not None:
                            cell_val = raw_val
                        r_cells[col_idx] = str(cell_val).strip()
                    elem.clear()
                    if r_cells:
                        max_c = max(r_cells.keys())
                        return [r_cells.get(i, "") for i in range(max_c + 1)]
    except Exception:
        return []
    return []


def _parse_full_xlsx_sheet(
    z: zipfile.ZipFile,
    target: str,
    shared_strings: list[str],
) -> tuple[list[str], list[list[Any]]]:
    """Parse all rows from a worksheet into (header, data_rows)."""
    if target not in z.namelist():
        return [], []
    headers: list[str] = []
    data_rows: list[list[Any]] = []

    try:
        with z.open(target) as stream:
            is_first = True
            for event, elem in ET.iterparse(stream, events=("end",)):
                if elem.tag == f"{{{NS_MAIN}}}row":
                    r_cells = {}
                    for c in elem.findall(f"{{{NS_MAIN}}}c"):
                        ref = c.attrib.get("r", "")
                        c_type = c.attrib.get("t", "")
                        _, col_idx = parse_cell_ref(ref)
                        v_elem = c.find(f"{{{NS_MAIN}}}v")
                        raw_val = v_elem.text if v_elem is not None else None
                        cell_val = None
                        if c_type == "s" and raw_val is not None:
                            s_idx = int(raw_val)
                            cell_val = shared_strings[s_idx] if s_idx < len(shared_strings) else ""
                        elif c_type == "inlineStr":
                            is_el = c.find(f"{{{NS_MAIN}}}is")
                            if is_el is not None:
                                t_el = is_el.find(f"{{{NS_MAIN}}}t")
                                if t_el is not None and t_el.text:
                                    cell_val = t_el.text
                        elif raw_val is not None:
                            cell_val = raw_val
                        r_cells[col_idx] = cell_val

                    if is_first:
                        if r_cells:
                            max_c = max(r_cells.keys())
                            headers = [str(r_cells.get(i, "") or "").strip() for i in range(max_c + 1)]
                        is_first = False
                    else:
                        if r_cells:
                            row_vals = [r_cells.get(i) for i in range(len(headers))]
                            data_rows.append(row_vals)
                    elem.clear()
    except Exception:
        return [], []

    return headers, data_rows


def _parse_rows_from_list(
    raw_data: list[list[Any]],
    headers: list[str],
    col_indices: dict[str, int],
) -> tuple[list[DegTableRow], list[AssetFinding]]:
    """Parse raw row arrays into structured DegTableRow instances with strict validation."""
    parsed_rows: list[DegTableRow] = []
    findings: list[AssetFinding] = []

    gene_idx = col_indices["gene_id"]
    symbol_idx = col_indices.get("gene_symbol")
    lfc_idx = col_indices["log2fc"]
    padj_idx = col_indices.get("padj")
    pval_idx = col_indices.get("pvalue")
    se_idx = col_indices.get("lfc_se")
    stat_idx = col_indices.get("stat")
    bm_idx = col_indices.get("base_mean")

    invalid_p_count = 0

    for row_idx, r in enumerate(raw_data, start=1):
        if not r or not any(v is not None and str(v).strip() != "" for v in r):
            continue

        def get_val(idx: int | None) -> Any:
            if idx is not None and idx < len(r):
                return r[idx]
            return None

        # 1. Gene identifier
        raw_gene = str(get_val(gene_idx) or "").strip()
        if not raw_gene:
            continue

        # 2. Gene symbol
        raw_symbol: str | None = None
        if symbol_idx is not None:
            s_val = get_val(symbol_idx)
            raw_symbol = str(s_val).strip() if s_val is not None and str(s_val).strip() != "" else None

        # 3. Log2 Fold Change
        lfc_val, lfc_unparsable = _parse_numeric(get_val(lfc_idx))

        # 4. P-value
        pval_val, _ = _parse_numeric(get_val(pval_idx))
        if pval_val is not None:
            if pval_val < -1e-7 or pval_val > 1.0 + 1e-7:
                invalid_p_count += 1
                pval_val = None

        # 5. Adjusted P-value / FDR
        padj_val, _ = _parse_numeric(get_val(padj_idx))
        if padj_val is not None:
            if padj_val < -1e-7 or padj_val > 1.0 + 1e-7:
                invalid_p_count += 1
                padj_val = None

        # 6. Optional fields - strictly None if not present or unparsable
        se_val, _ = _parse_numeric(get_val(se_idx))
        stat_val, _ = _parse_numeric(get_val(stat_idx))
        bm_val, _ = _parse_numeric(get_val(bm_idx))

        # Extra fields preserved for provenance
        extra = {}
        for c_idx, val in enumerate(r):
            if c_idx not in {gene_idx, symbol_idx, lfc_idx, padj_idx, pval_idx, se_idx, stat_idx, bm_idx}:
                if c_idx < len(headers):
                    extra[headers[c_idx]] = val

        parsed_rows.append(
            DegTableRow(
                source_index=row_idx,
                raw_gene_id=raw_gene,
                raw_gene_symbol=raw_symbol,
                log2_fold_change=lfc_val,
                p_value=pval_val,
                adjusted_p_value=padj_val,
                lfc_standard_error=se_val,
                wald_statistic=stat_val,
                base_mean=bm_val,
                extra_fields=extra,
            )
        )

    if invalid_p_count > 0:
        findings.append(
            AssetFinding(
                "WARNING",
                "invalid_significance_values",
                "deg_table",
                f"Encountered {invalid_p_count} p-value or padj values outside the valid [0, 1] range; set to None.",
            )
        )

    return parsed_rows, findings
