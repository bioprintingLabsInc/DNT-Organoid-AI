"""Safe, bounded, and deterministic Excel (.xlsx) workbook inspection and parsing.

Implements pure-Python parsing via standard library zipfile and xml.etree.ElementTree
with zero external dependencies. Enforces:
1. Rejection of corrupt or password-protected workbooks.
2. Safe sheet enumeration.
3. Candidate expression sheet identification without silent multi-sheet selection.
4. Formula rejection unless cached numeric value is verified.
5. Exact integer representation preservation.
6. Provenance capture of workbook path, selected sheet, and formula presence.
"""

from collections import Counter
import re
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

from src.ingestion.models import AssetFinding
from src.qc.matrix import CountMatrix
from src.qc.result import Finding, Severity

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def col2num(col_str: str) -> int:
    """Convert spreadsheet column letters (A, B, ..., Z, AA, ...) to 0-based index."""
    num = 0
    for c in col_str:
        num = num * 26 + (ord(c.upper()) - ord("A")) + 1
    return num - 1


def parse_cell_ref(ref: str) -> tuple[int, int]:
    """Parse cell reference (e.g. 'B12') into (row_idx, col_idx), 0-based."""
    m = re.match(r"^([A-Za-z]+)(\d+)$", ref)
    if not m:
        return 0, 0
    col_str, row_str = m.groups()
    return int(row_str) - 1, col2num(col_str)


def inspect_xlsx_workbook(path: Path) -> tuple[bool, str | None, list[dict[str, str]]]:
    """Inspect workbook integrity, password protection, and enumerate sheets.

    Returns: (is_valid, error_message, sheets_info)
    """
    if not path.is_file():
        return False, f"File not found: {path}", []
    if not zipfile.is_zipfile(path):
        return False, "File is corrupt or not a valid XLSX archive", []

    try:
        with zipfile.ZipFile(path, "r") as z:
            # Check for password protection / encryption in any entry
            for info in z.infolist():
                if info.flag_bits & 0x1:
                    return False, "Workbook is password-protected or encrypted", []

            namelist = z.namelist()
            if "xl/workbook.xml" not in namelist:
                return False, "Malformed XLSX: missing xl/workbook.xml", []

            # Parse rels to map rId -> target sheet xml
            rels_map = {}
            if "xl/_rels/workbook.xml.rels" in namelist:
                rels_root = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
                for elem in rels_root:
                    r_id = elem.attrib.get("Id")
                    target = elem.attrib.get("Target", "")
                    if target and not target.startswith("xl/"):
                        if not target.startswith("/"):
                            target = "xl/" + target
                    rels_map[r_id] = target

            wb_root = ET.fromstring(z.read("xl/workbook.xml"))
            sheets_elem = wb_root.find(f"{{{NS_MAIN}}}sheets")
            if sheets_elem is None:
                return False, "Malformed XLSX: no sheets element found in workbook.xml", []

            sheets = []
            for s in sheets_elem:
                name = s.attrib.get("name", "")
                r_id = s.attrib.get(f"{{{NS_REL}}}id", "")
                target = rels_map.get(r_id, "")
                sheets.append({"name": name, "rId": r_id, "target": target})

            return True, None, sheets
    except Exception as e:
        return False, f"Failed to inspect XLSX workbook: {e}", []


def read_shared_strings(z: zipfile.ZipFile) -> list[str]:
    """Read xl/sharedStrings.xml if present."""
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    shared: list[str] = []
    try:
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall(f"{{{NS_MAIN}}}si"):
            t = si.find(f"{{{NS_MAIN}}}t")
            if t is not None and t.text:
                shared.append(t.text)
            else:
                # Handle formatted runs <r><t>...</t></r>
                runs = si.findall(f"{{{NS_MAIN}}}r")
                if runs:
                    text_parts = [r.find(f"{{{NS_MAIN}}}t").text for r in runs if r.find(f"{{{NS_MAIN}}}t") is not None and r.find(f"{{{NS_MAIN}}}t").text]
                    shared.append("".join(text_parts))
                else:
                    shared.append("")
    except Exception:
        pass
    return shared


def inspect_sheet_candidate(
    z: zipfile.ZipFile,
    target: str,
    shared_strings: list[str],
    max_rows: int = 51,
) -> dict[str, Any]:
    """Quickly inspect up to max_rows of a worksheet to evaluate if it is a count matrix."""
    if target not in z.namelist():
        return {"is_candidate": False, "reason": "Target sheet XML not found"}

    rows_data: dict[int, dict[int, Any]] = {}
    formulas_present = False
    has_unevaluated_formula = False

    try:
        with z.open(target) as stream:
            for event, elem in ET.iterparse(stream, events=("end",)):
                if elem.tag == f"{{{NS_MAIN}}}row":
                    r_attr = elem.attrib.get("r")
                    row_idx = int(r_attr) - 1 if r_attr else len(rows_data)
                    row_cells = {}
                    for c in elem.findall(f"{{{NS_MAIN}}}c"):
                        ref = c.attrib.get("r", "")
                        c_type = c.attrib.get("t", "")
                        _, col_idx = parse_cell_ref(ref)

                        f_elem = c.find(f"{{{NS_MAIN}}}f")
                        v_elem = c.find(f"{{{NS_MAIN}}}v")

                        if f_elem is not None:
                            formulas_present = True
                            if v_elem is None or not v_elem.text:
                                has_unevaluated_formula = True

                        raw_val = v_elem.text if v_elem is not None else None
                        cell_val = None
                        if c_type == "inlineStr":
                            is_elem = c.find(f"{{{NS_MAIN}}}is")
                            if is_elem is not None:
                                t_elem = is_elem.find(f"{{{NS_MAIN}}}t")
                                cell_val = t_elem.text if t_elem is not None and t_elem.text else ""
                        elif c_type == "s" and raw_val is not None:
                            s_idx = int(raw_val)
                            cell_val = shared_strings[s_idx] if s_idx < len(shared_strings) else ""
                        elif c_type == "b" and raw_val is not None:
                            cell_val = raw_val == "1"
                        elif raw_val is not None:
                            cell_val = raw_val
                        row_cells[col_idx] = cell_val

                    rows_data[row_idx] = row_cells
                    elem.clear()
                    if len(rows_data) >= max_rows:
                        break
    except Exception as e:
        return {"is_candidate": False, "reason": f"XML parse error: {e}"}

    if len(rows_data) < 2:
        return {"is_candidate": False, "reason": "Fewer than 2 rows in sheet"}

    sorted_row_indices = sorted(rows_data.keys())
    header_row = rows_data[sorted_row_indices[0]]
    all_col_indices = set()
    for r in sorted_row_indices:
        all_col_indices.update(rows_data[r].keys())
    if not all_col_indices:
        return {"is_candidate": False, "reason": "No cells found"}

    min_col = min(all_col_indices)
    max_col = max(all_col_indices)
    col_count = max_col - min_col + 1
    if col_count < 2:
        return {"is_candidate": False, "reason": "Fewer than 2 columns"}

    # Check data rows for non-negative integers in expression columns
    data_values: list[float] = []
    has_strings_in_data = False
    for r_idx in sorted_row_indices[1:]:
        row_dict = rows_data[r_idx]
        for c_idx in range(min_col + 1, max_col + 1):
            val = row_dict.get(c_idx)
            if val is None or str(val).strip() == "":
                continue
            try:
                num = float(str(val).strip())
                data_values.append(num)
            except ValueError:
                has_strings_in_data = True
                break
        if has_strings_in_data:
            break

    if has_strings_in_data or not data_values:
        return {"is_candidate": False, "reason": "Data columns are not numeric"}

    is_integer_like = all(v >= 0 and v.is_integer() for v in data_values)
    return {
        "is_candidate": is_integer_like,
        "is_integer_like": is_integer_like,
        "data_count": len(data_values),
        "formulas_present": formulas_present,
        "has_unevaluated_formula": has_unevaluated_formula,
        "min_col": min_col,
        "max_col": max_col,
    }


def classify_xlsx_workbook(path: Path) -> tuple[str | None, list[str], list[AssetFinding], dict[str, Any]]:
    """Determine whether workbook contains exactly one valid count matrix sheet.

    Returns: (selected_sheet_name, all_sheets, findings, provenance)
    """
    valid, err, sheets = inspect_xlsx_workbook(path)
    if not valid:
        finding = AssetFinding("ERROR", "corrupt_workbook", "source_reference", f"Invalid XLSX workbook: {err}")
        return None, [], [finding], {"status": "CORRUPT_OR_ENCRYPTED", "error": err}

    all_sheet_names = [s["name"] for s in sheets]
    if not sheets:
        finding = AssetFinding("ERROR", "empty_workbook", "source_reference", "Workbook contains no worksheets.")
        return None, [], [finding], {"status": "EMPTY_WORKBOOK"}

    candidates = []
    sheet_details = {}

    with zipfile.ZipFile(path, "r") as z:
        shared = read_shared_strings(z)
        for s in sheets:
            target = s["target"]
            c_info = inspect_sheet_candidate(z, target, shared)
            sheet_details[s["name"]] = c_info
            if c_info["is_candidate"] and not c_info.get("has_unevaluated_formula", False):
                candidates.append(s["name"])

    provenance = {
        "workbook_filename": path.name,
        "all_sheets": all_sheet_names,
        "candidate_sheets": candidates,
    }

    if len(candidates) == 1:
        selected = candidates[0]
        provenance["selected_sheet"] = selected
        provenance["formulas_present"] = sheet_details[selected].get("formulas_present", False)
        return selected, all_sheet_names, [], provenance
    elif len(candidates) > 1:
        finding = AssetFinding(
            "WARNING",
            "ambiguous_sheets",
            "source_reference",
            f"Multiple candidate expression sheets detected ({candidates}); manual review is required.",
        )
        provenance["status"] = "AMBIGUOUS_CANDIDATE_SHEETS"
        return None, all_sheet_names, [finding], provenance
    else:
        finding = AssetFinding(
            "WARNING",
            "no_valid_count_sheet",
            "source_reference",
            "No sheet containing a valid non-negative integer count matrix was identified.",
        )
        provenance["status"] = "NO_VALID_COUNT_SHEET"
        return None, all_sheet_names, [finding], provenance


def read_xlsx_count_matrix(
    path: str | Path,
    rules: dict | None = None,
    sheet_name: str | None = None,
) -> tuple[CountMatrix | None, list[Finding], dict[str, Any]]:
    """Parse a valid count matrix sheet from an Excel workbook into a CountMatrix.

    Canonical output is equivalent to CountMatrix with full Stage 6A compatibility.
    """
    path_p = Path(path)
    valid, err, sheets = inspect_xlsx_workbook(path_p)
    if not valid:
        return None, [Finding(Severity.ERROR, "corrupt_workbook", "dataset", None, "matrix", f"Cannot open XLSX: {err}")], {}

    target_sheet = sheet_name
    all_sheet_names = [s["name"] for s in sheets]
    provenance: dict[str, Any] = {
        "workbook_filename": path_p.name,
        "all_sheets": all_sheet_names,
    }

    if not target_sheet:
        selected, _, findings_asset, prov = classify_xlsx_workbook(path_p)
        provenance.update(prov)
        if selected is None:
            qc_findings = [
                Finding(
                    Severity.ERROR if f.severity == "ERROR" else Severity.REVIEW,
                    f.rule,
                    "dataset",
                    None,
                    "sheet_selection",
                    f.message,
                )
                for f in findings_asset
            ]
            return None, qc_findings, provenance
        target_sheet = selected

    provenance["selected_sheet"] = target_sheet
    target_info = next((s for s in sheets if s["name"] == target_sheet), None)
    if not target_info:
        return None, [Finding(Severity.ERROR, "sheet_not_found", "dataset", None, "sheet", f"Sheet '{target_sheet}' not found.")], provenance

    findings: list[Finding] = []
    rows_data: dict[int, dict[int, Any]] = {}
    formulas_present = False
    has_unevaluated_formula = False

    try:
        with zipfile.ZipFile(path_p, "r") as z:
            shared = read_shared_strings(z)
            with z.open(target_info["target"]) as stream:
                for event, elem in ET.iterparse(stream, events=("end",)):
                    if elem.tag == f"{{{NS_MAIN}}}row":
                        r_attr = elem.attrib.get("r")
                        row_idx = int(r_attr) - 1 if r_attr else len(rows_data)
                        row_cells = {}
                        for c in elem.findall(f"{{{NS_MAIN}}}c"):
                            ref = c.attrib.get("r", "")
                            c_type = c.attrib.get("t", "")
                            _, col_idx = parse_cell_ref(ref)

                            f_elem = c.find(f"{{{NS_MAIN}}}f")
                            v_elem = c.find(f"{{{NS_MAIN}}}v")

                            if f_elem is not None:
                                formulas_present = True
                                if v_elem is None or not v_elem.text:
                                    has_unevaluated_formula = True

                            raw_val = v_elem.text if v_elem is not None else None
                            cell_val = None
                            if c_type == "inlineStr":
                                is_elem = c.find(f"{{{NS_MAIN}}}is")
                                if is_elem is not None:
                                    t_elem = is_elem.find(f"{{{NS_MAIN}}}t")
                                    cell_val = t_elem.text if t_elem is not None and t_elem.text else ""
                            elif c_type == "s" and raw_val is not None:
                                s_idx = int(raw_val)
                                cell_val = shared[s_idx] if s_idx < len(shared) else ""
                            elif c_type == "b" and raw_val is not None:
                                cell_val = raw_val == "1"
                            elif raw_val is not None:
                                cell_val = raw_val
                            row_cells[col_idx] = cell_val

                        rows_data[row_idx] = row_cells
                        elem.clear()
    except Exception as e:
        return None, [Finding(Severity.ERROR, "malformed_matrix", "dataset", None, "matrix", f"Error reading sheet XML: {e}")], provenance

    provenance["formulas_present"] = formulas_present
    if has_unevaluated_formula:
        findings.append(Finding(Severity.ERROR, "unevaluated_formulas_forbidden", "dataset", None, "cells", "Sheet contains formulas without cached numeric values."))
        return None, findings, provenance

    if len(rows_data) < 2:
        return None, [Finding(Severity.ERROR, "malformed_matrix", "dataset", None, "matrix", "Sheet contains fewer than 2 rows.")], provenance

    sorted_row_indices = sorted(rows_data.keys())
    first_row = rows_data[sorted_row_indices[0]]

    # Detect bounding columns
    all_cols = set()
    for r in rows_data.values():
        all_cols.update(r.keys())
    if not all_cols:
        return None, [Finding(Severity.ERROR, "malformed_matrix", "dataset", None, "matrix", "Sheet has no cells.")], provenance

    min_col = min(all_cols)
    max_col = max(all_cols)
    if (max_col - min_col) < 1:
        return None, [Finding(Severity.ERROR, "sample_columns_required", "dataset", None, "header", "At least one sample column is required.")], provenance

    # Check sample headers in row 1
    sample_col_indices = list(range(min_col + 1, max_col + 1))
    sample_ids = []
    for c_idx in sample_col_indices:
        h_val = first_row.get(c_idx)
        s_name = str(h_val).strip() if h_val is not None else f"Sample_{c_idx}"
        sample_ids.append(s_name)

    # Detect gene column header
    gene_header_val = first_row.get(min_col)
    gene_col_header = str(gene_header_val).strip() if gene_header_val is not None and str(gene_header_val).strip() else "gene_id"
    provenance["gene_column_header"] = gene_col_header

    # Duplicate samples check (O(N))
    sample_counts = Counter(sample_ids)
    for s in sorted([x for x, c in sample_counts.items() if c > 1]):
        findings.append(Finding(Severity.ERROR, "unique_sample_identifiers_required", "sample", s, "header", "Duplicate sample identifier."))

    genes: list[str] = []
    values: list[list[int]] = [[] for _ in sample_ids]

    for line_idx, r_idx in enumerate(sorted_row_indices[1:], start=2):
        row_dict = rows_data[r_idx]
        g_val = row_dict.get(min_col)
        g_str = str(g_val).strip() if g_val is not None else ""
        if not g_str:
            findings.append(Finding(Severity.ERROR, "empty_gene_identifiers", "gene", None, f"row[{line_idx}].gene_id", "Empty gene identifier."))
            continue
        genes.append(g_str)

        for j, c_idx in enumerate(sample_col_indices):
            cell_val = row_dict.get(c_idx)
            if cell_val is None or str(cell_val).strip() == "":
                findings.append(Finding(Severity.ERROR, "missing_counts_forbidden", "sample", sample_ids[j], f"row[{line_idx}].{sample_ids[j]}", "Missing count value."))
                continue
            try:
                num = float(str(cell_val).strip())
            except ValueError:
                findings.append(Finding(Severity.ERROR, "numeric_counts_required", "sample", sample_ids[j], f"row[{line_idx}].{sample_ids[j]}", "Count must be numeric."))
                continue
            if num < 0:
                findings.append(Finding(Severity.ERROR, "non_negative_counts_required", "sample", sample_ids[j], f"row[{line_idx}].{sample_ids[j]}", "Count must be non-negative."))
            if not num.is_integer():
                findings.append(Finding(Severity.ERROR, "integer_like_counts_required", "sample", sample_ids[j], f"row[{line_idx}].{sample_ids[j]}", "Count must be integer-like."))
            values[j].append(int(num))

    if any(f.severity == Severity.ERROR for f in findings):
        return None, findings, provenance

    # Duplicate genes check (O(N) Counter)
    gene_counts = Counter(genes)
    for gene in sorted([x for x, c in gene_counts.items() if c > 1]):
        findings.append(Finding(Severity.REVIEW, "duplicate_gene_identifiers", "gene", gene, "gene_id", "Duplicate gene identifier; no aggregation performed."))

    # Zero only genes
    for i, gene in enumerate(genes):
        if values and all(col[i] == 0 for col in values):
            findings.append(Finding(Severity.INFO, "zero_only_genes", "gene", gene, f"gene[{i}]", "Gene has zero counts in every sample."))

    # Zero only samples
    for sample, col in zip(sample_ids, values):
        if col and not any(col):
            findings.append(Finding(Severity.REVIEW, "zero_only_samples", "sample", sample, f"sample[{sample}]", "Sample has zero counts for every gene."))

    matrix = CountMatrix(tuple(genes), tuple(sample_ids), tuple(tuple(v) for v in values))
    return matrix, findings, provenance
