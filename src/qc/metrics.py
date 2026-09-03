"""Deterministic descriptive count metrics."""
from statistics import median

def quantile(values, q):
    ordered=sorted(values); position=(len(ordered)-1)*q; lo=int(position); hi=min(lo+1,len(ordered)-1); fraction=position-lo
    return ordered[lo]*(1-fraction)+ordered[hi]*fraction

def basic_metrics(values):
    total=sum(values); zero=sum(v==0 for v in values)
    return {"library_size":total,"detected_genes":len(values)-zero,"zero_count_genes":zero,"zero_fraction":zero/len(values) if values else 0.0,"minimum":min(values,default=0),"median":median(values) if values else 0.0,"maximum":max(values,default=0),"quantiles":(quantile(values,.25),quantile(values,.5),quantile(values,.75)) if values else (0.0,0.0,0.0)}

def ratio(value, reference): return value/reference if reference else None
