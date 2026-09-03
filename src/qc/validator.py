"""Step 6A bulk raw-count QC orchestration."""
from dataclasses import replace
from pathlib import Path
from statistics import median
from src.ingestion.registry import reported_value
from .config import load_rules,rules_checksum
from .matrix import read_matrix
from .metadata import reconcile
from .metrics import basic_metrics,ratio
from .replicates import condition_results
from .result import DatasetResult,Finding,SampleResult,Severity,ordered,status_for

CALCULATION_VERSION="bulk_raw_count_qc_v1"

def validate_bulk_counts(path:str|Path,metadata:dict,asset_id:str,assay_id:str,rules_path=None):
    rules=load_rules(rules_path) if rules_path else load_rules(); matrix,findings=read_matrix(path,rules)
    if matrix is None: return DatasetResult(status_for(findings),ordered(findings),0,0,asset_id,assay_id,(),(),rules["rules_version"],rules_checksum(rules),CALCULATION_VERSION,(),())
    contexts,metadata_findings=reconcile(metadata,matrix.sample_ids,asset_id,assay_id); findings+=metadata_findings
    raw=[basic_metrics(col) for col in matrix.columns]; dataset_median=median([m["library_size"] for m in raw]) if raw else 0
    condition_medians={}
    for sid,m in zip(matrix.sample_ids,raw):
        condition=contexts.get(sid,{}).get("condition_id"); condition_medians.setdefault(condition,[]).append(m["library_size"])
    condition_medians={k:median(v) for k,v in condition_medians.items()}
    samples=[]
    for sid,m in zip(matrix.sample_ids,raw):
        context=contexts.get(sid,{}); samples.append(SampleResult(sid,context.get("condition_id"),context.get("biological_replicate_id"),context.get("technical_replicate_id"),m["library_size"],m["detected_genes"],m["zero_count_genes"],m["zero_fraction"],m["minimum"],m["median"],m["maximum"],m["quantiles"],ratio(m["library_size"],dataset_median),ratio(m["library_size"],condition_medians.get(context.get("condition_id"),0)),()))
    from .rules import evaluate_thresholds
    for i,sample in enumerate(samples):
        sf=evaluate_thresholds(sample,rules["warning_thresholds"]); findings+=sf; samples[i]=replace(sample,findings=ordered(sf))
    conditions=condition_results(tuple(samples),matrix.columns,rules["replicate_correlation"]["implementation"])
    findings += [finding for condition in conditions for finding in condition.findings]
    studies=tuple(sorted({e.get("study_id") for e in metadata.get("experiments",[]) if e.get("study_id")})); experiments=tuple(sorted({c.get("experiment_id") for c in contexts.values() if c.get("experiment_id")}))
    findings=ordered(findings)
    return DatasetResult(status_for(findings),findings,len(matrix.gene_ids),len(matrix.sample_ids),asset_id,assay_id,studies,experiments,rules["rules_version"],rules_checksum(rules),CALCULATION_VERSION,tuple(samples),conditions)
