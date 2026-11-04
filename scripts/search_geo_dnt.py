#!/usr/bin/env python3
"""Reproducible NCBI GEO Series discovery. Network access required.

Usage:
    python scripts/search_geo_dnt.py --output-dir outputs/geo_search/RUN_ID

Records exact ESearch queries, NCBI counts, retrieved UID/accessions, API
responses metadata, UTC timestamps, and deduplicated accession membership.
No historical collection claims are inferred from a new run.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import pathlib
import platform
import sys
import time
import urllib.parse
import urllib.request

API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def get(endpoint, params):
    url = API + endpoint + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "DNT-Organoid-AI-GEO-discovery/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                return json.load(response)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")

def or_group(terms):
    # Quoted phrases in All Fields; NCBI interprets phrases according to indexing.
    return "(" + " OR ".join('"' + t.replace('"', '') + '"[All Fields]' for t in terms) + ")"

def build_query(cfg, q):
    groups = cfg["groups"]
    left = [s for g in q["join"] for s in groups[g]]
    right = groups[q["with"]]
    return or_group(left) + " AND " + or_group(right) + " AND " + " AND ".join(cfg["filters"])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True)
    p.add_argument("--config", default="config/dnt_geo_search_keywords_v1.json")
    p.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
    p.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    args = p.parse_args()
    if not args.email:
        p.error("Pass --email or set NCBI_EMAIL to identify requests to NCBI.")
    output = pathlib.Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=False)
    cfgpath = pathlib.Path(args.config)
    cfg = json.loads(cfgpath.read_text(encoding="utf-8"))
    run = {"protocol":cfg["protocol_id"],"run_started_utc":now(),
           "config_sha256":hashlib.sha256(cfgpath.read_bytes()).hexdigest(),
           "python":sys.version,"platform":platform.platform(),"queries":[]}
    memberships = {}
    for q in cfg["queries"]:
        query = build_query(cfg,q)
        start = now()
        base = {"db":"gds","term":query,"retmode":"json",
                "tool":"DNT-Organoid-AI-GEO-discovery","email":args.email}
        if args.api_key:
            base["api_key"] = args.api_key
        first = get("esearch.fcgi",{**base,"retmax":0})
        total = int(first["esearchresult"]["count"])
        uids = []
        for offset in range(0,total,500):
            page = get("esearch.fcgi",{**base,"retstart":offset,"retmax":500})
            uids.extend(page["esearchresult"]["idlist"])
            time.sleep(0.4)
        if len(set(uids)) != total:
            raise RuntimeError(f"{q['id']}: ESearch count/UID mismatch; rerun required")
        accessions = []
        for offset in range(0,len(uids),100):
            batch = uids[offset:offset+100]
            summary_params = {"db":"gds","id":",".join(batch),
                        "retmode":"json","tool":"DNT-Organoid-AI-GEO-discovery",
                        "email":args.email}
            if args.api_key:
                summary_params["api_key"] = args.api_key
            summary = get("esummary.fcgi", summary_params)
            for uid in batch:
                record = summary["result"].get(str(uid),{})
                acc = str(record.get("accession",""))
                if not acc.startswith("GSE"):
                    raise RuntimeError(f"Missing/non-GSE accession for UID {uid}: {acc}")
                accessions.append(acc)
                memberships.setdefault(acc,set()).add(q["id"])
            time.sleep(0.4)
        with (output/(q["id"]+"_accessions.csv")).open("w",newline="",encoding="utf-8") as f:
            writer=csv.writer(f);writer.writerow(["geo_accession"])
            writer.writerows([[a] for a in accessions])
        run["queries"].append({"query_id":q["id"],"exact_query":query,
                                "start_utc":start,"end_utc":now(),
                                "esearch_reported_count":total,
                                "retrieved_gse_accessions":len(accessions)})
    with (output/"all_unique_geo_series.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f);writer.writerow(["geo_accession","matched_query_ids"])
        for acc in sorted(memberships):
            writer.writerow([acc,";".join(sorted(memberships[acc]))])
    run["unique_accessions"]=len(memberships)
    run["run_finished_utc"]=now()
    (output/"run_manifest.json").write_text(json.dumps(run,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"run_finished_utc":run["run_finished_utc"],
                      "unique_accessions":run["unique_accessions"],
                      "query_counts":[(q["query_id"],q["esearch_reported_count"]) for q in run["queries"]]},indent=2))

if __name__ == "__main__":
    main()
