"""Canonical metadata reconciliation; never infers biological roles."""
from src.ingestion.registry import reported_value
from .result import Finding,Severity

def reconcile(metadata,matrix_samples,asset_id,assay_id):
    findings=[]; samples=[s for s in metadata.get("samples",[]) if isinstance(s,dict)]; index={}
    for s in samples: index.setdefault(s.get("sample_id"),[]).append(s)
    assay=next((a for a in metadata.get("sequencing_assays",[]) if a.get("assay_id")==assay_id),None)
    asset=next((a for a in metadata.get("input_data_assets",[]) if a.get("asset_id")==asset_id),None)
    if not assay: findings.append(Finding(Severity.ERROR,"assay_missing","assay",assay_id,"assay_id","Assay not found.")); expected=set()
    else: expected=set(assay.get("sample_ids",[]))
    if not asset or assay_id not in asset.get("assay_ids",[]): findings.append(Finding(Severity.ERROR,"asset_assay_mismatch","asset",asset_id,"assay_ids","Asset does not belong to assay."))
    for sid in matrix_samples:
        if sid not in index: findings.append(Finding(Severity.REVIEW,"unexpected_matrix_samples","sample",sid,"matrix.header","Matrix sample is absent from canonical metadata."))
        elif len(index[sid])!=1: findings.append(Finding(Severity.ERROR,"ambiguous_sample_mapping","sample",sid,"samples","Sample does not map uniquely."))
        elif sid not in expected: findings.append(Finding(Severity.ERROR,"assay_sample_mismatch","sample",sid,"assay.sample_ids","Mapped sample is not assigned to assay."))
    for sid in sorted(expected-set(matrix_samples)): findings.append(Finding(Severity.REVIEW,"missing_expected_samples","sample",sid,"matrix.header","Expected assay sample is absent from matrix."))
    contexts={s.get("sample_id"): {"condition_id":s.get("condition_id"),"experiment_id":s.get("experiment_id"),"biological_replicate_id":reported_value(s.get("biological_replicate_id")),"technical_replicate_id":reported_value(s.get("technical_replicate_id"))} for s in samples if s.get("sample_id") in matrix_samples}
    return contexts,findings
