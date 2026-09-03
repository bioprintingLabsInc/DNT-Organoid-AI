"""Generic opt-in numerical threshold evaluation."""
from .result import Finding,Severity

def evaluate_thresholds(sample,thresholds):
    findings=[]; mapping={"minimum_library_size":("library_size",lambda a,b:a<b),"minimum_detected_genes":("detected_genes",lambda a,b:a<b),"maximum_zero_gene_fraction":("zero_fraction",lambda a,b:a>b),"relative_library_size_lower_ratio":("dataset_library_ratio",lambda a,b:a is not None and a<b),"relative_library_size_upper_ratio":("dataset_library_ratio",lambda a,b:a is not None and a>b)}
    for key,(field,trigger) in mapping.items():
        threshold=thresholds.get(key)
        if threshold is not None and trigger(getattr(sample,field),threshold): findings.append(Finding(Severity.WARNING,key,"sample",sample.sample_id,field,f"Configured threshold {threshold} triggered."))
    return findings
