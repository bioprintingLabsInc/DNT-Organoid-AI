"""Non-mutating bulk count-matrix parsing and structural QC."""
import csv
from dataclasses import dataclass
from pathlib import Path
from .result import Finding, Severity

@dataclass(frozen=True, slots=True)
class CountMatrix:
    gene_ids: tuple[str,...]; sample_ids: tuple[str,...]; columns: tuple[tuple[int,...],...]

def read_matrix(path: str|Path, rules: dict) -> tuple[CountMatrix|None,list[Finding]]:
    findings=[]; path=Path(path)
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.reader(stream, delimiter="," if path.suffix.lower()==".csv" else "\t"))
    except (OSError,UnicodeError,csv.Error) as e: return None,[Finding(Severity.ERROR,"malformed_matrix","dataset",None,"matrix",str(e))]
    if not rows or not any(rows): return None,[Finding(Severity.ERROR,"malformed_matrix","dataset",None,"matrix","Matrix is empty.")]
    header=rows[0]
    if not header or not header[0].strip(): findings.append(Finding(Severity.ERROR,"gene_identifier_column_required","dataset",None,"header[0]","Gene identifier column is missing."))
    samples=tuple(header[1:])
    if not samples: findings.append(Finding(Severity.ERROR,"sample_columns_required","dataset",None,"header","At least one sample column is required."))
    for s in sorted({x for x in samples if samples.count(x)>1}): findings.append(Finding(Severity.ERROR,"unique_sample_identifiers_required","sample",s,"header","Duplicate sample identifier."))
    genes=[]; values=[[] for _ in samples]
    for i,row in enumerate(rows[1:],1):
        if not row or not any(cell.strip() for cell in row): findings.append(Finding(Severity.ERROR,"empty_row","gene",None,f"row[{i}]","Empty row.")); continue
        if len(row)!=len(header): findings.append(Finding(Severity.ERROR,"malformed_row","gene",row[0] if row else None,f"row[{i}]","Row width differs from header.")); continue
        gene=row[0].strip(); genes.append(gene)
        if not gene: findings.append(Finding(Severity.ERROR,"empty_gene_identifiers","gene",None,f"row[{i}].gene_id","Empty gene identifier."))
        for j,text in enumerate(row[1:]):
            if text.strip()=="": findings.append(Finding(Severity.ERROR,"missing_counts_forbidden","sample",samples[j],f"row[{i}].{samples[j]}","Missing count value.")); continue
            try: number=float(text)
            except ValueError: findings.append(Finding(Severity.ERROR,"numeric_counts_required","sample",samples[j],f"row[{i}].{samples[j]}","Count must be numeric.")); continue
            if number<0: findings.append(Finding(Severity.ERROR,"non_negative_counts_required","sample",samples[j],f"row[{i}].{samples[j]}","Count must be non-negative."))
            if not number.is_integer(): findings.append(Finding(Severity.ERROR,"integer_like_counts_required","sample",samples[j],f"row[{i}].{samples[j]}","Count must be integer-like."))
            values[j].append(int(number))
    if any(f.severity==Severity.ERROR for f in findings): return None,findings
    for gene in sorted({x for x in genes if genes.count(x)>1}): findings.append(Finding(Severity.REVIEW,"duplicate_gene_identifiers","gene",gene,"gene_id","Duplicate gene identifier; no aggregation performed."))
    for i,gene in enumerate(genes):
        if values and all(col[i]==0 for col in values): findings.append(Finding(Severity.INFO,"zero_only_genes","gene",gene,f"gene[{i}]","Gene has zero counts in every sample."))
    for sample,col in zip(samples,values):
        if col and not any(col): findings.append(Finding(Severity.REVIEW,"zero_only_samples","sample",sample,f"sample[{sample}]","Sample has zero counts for every gene."))
    return CountMatrix(tuple(genes),samples,tuple(tuple(v) for v in values)),findings
