"""Condition-level biological replicate Spearman summaries."""
from itertools import combinations
from math import sqrt
from statistics import median
from .result import ConditionResult,Finding,PairCorrelation,Severity

def ranks(values):
    out=[0.0]*len(values); ordered=sorted(enumerate(values),key=lambda x:x[1]); i=0
    while i<len(ordered):
        j=i
        while j+1<len(ordered) and ordered[j+1][1]==ordered[i][1]: j+=1
        rank=(i+j+2)/2
        for k in range(i,j+1): out[ordered[k][0]]=rank
        i=j+1
    return out

def spearman(a,b):
    ra,rb=ranks(a),ranks(b); ma,mb=sum(ra)/len(ra),sum(rb)/len(rb); da=[x-ma for x in ra]; db=[x-mb for x in rb]
    denominator=sqrt(sum(x*x for x in da)*sum(x*x for x in db))
    return None if denominator==0 else sum(x*y for x,y in zip(da,db))/denominator

def condition_results(sample_results,columns,implementation):
    groups={}
    for result,col in zip(sample_results,columns): groups.setdefault(result.condition_id,[]).append((result,col))
    output=[]
    for condition,items in sorted(groups.items(),key=lambda x:x[0] or ""):
        by_rep={}; technical={}
        for result,col in items:
            by_rep.setdefault(result.biological_replicate_id,[]).append((result,col))
            if result.technical_replicate_id is not None: technical.setdefault(result.biological_replicate_id,[]).append(result.technical_replicate_id)
        pairs=[]; reasons=[]; findings=[]
        for (ra,aa),(rb,bb) in combinations(sorted(by_rep.items(),key=lambda x:x[0] or ""),2):
            if len(aa)!=1 or len(bb)!=1:
                value=None; reason="technical_replicates_require_unapproved_aggregation"; reasons.append(reason)
                findings.append(Finding(Severity.REVIEW,reason,"condition",condition or None,"technical_replicates","Technical replicates make a biological-replicate profile ambiguous; no aggregation or selection was performed."))
            else: value=spearman(aa[0][1],bb[0][1]); reason=None if value is not None else "constant_or_all_zero_vector"; reasons += [reason] if reason else []
            pairs.append(PairCorrelation(aa[0][0].sample_id,bb[0][0].sample_id,ra or "",rb or "",value,reason))
        evaluable=[p.value for p in pairs if p.value is not None]
        summary=(min(evaluable),median(evaluable),max(evaluable)) if evaluable else None
        if len(by_rep)<2: reasons.append("fewer_than_two_biological_replicates")
        output.append(ConditionResult(condition or "",len(items),len(by_rep),tuple((k or "",tuple(sorted(v))) for k,v in sorted(technical.items(),key=lambda x:x[0] or "")),"spearman",implementation,tuple(pairs),len(evaluable),summary,tuple(sorted(set(reasons))),tuple(findings)))
    return tuple(output)
