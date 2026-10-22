"""Exhaustive live re-verification of the 12 approved PRJNA802072 pilot runs against ENA/SRA API."""

import csv
import json
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPO_ROOT / "config" / "quantification_bridge" / "pilot_manifest.json"

print(f"Loading governed manifest: {MANIFEST_PATH}")
with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
    manifest = json.load(f)

url = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?"
    "accession=PRJNA802072&result=read_run&"
    "fields=study_accession,sample_accession,experiment_accession,run_accession,"
    "sample_alias,sample_title,fastq_ftp,fastq_md5,read_count,base_count,library_layout&"
    "format=tsv&download=true"
)

print(f"Querying live ENA API: {url}...")
req = urllib.request.Request(url, headers={"User-Agent": "DNT-Organoid-AI/2.0"})
with urllib.request.urlopen(req) as resp:
    lines = resp.read().decode("utf-8").strip().split("\n")

reader = csv.DictReader(lines, delimiter="\t")
ena_map = {r["run_accession"]: r for r in reader}
print(f"Total runs returned by ENA for PRJNA802072: {len(ena_map)}")

approved_runs = manifest["samples"]
assert len(approved_runs) == 12, f"Expected exactly 12 runs, got {len(approved_runs)}"

print("\n" + "=" * 120)
print(f"{'Run ID':<13} {'Sample Alias':<13} {'Group':<10} {'Dose':<8} {'Day':<4} {'Layout':<8} {'BioSample':<14} {'Experiment':<13} {'Read Pairs':<12} {'ENA MD5 Match'}")
print("=" * 120)

all_reconciled = True
for s in approved_runs:
    rid = s["run_accession"]
    if rid not in ena_map:
        print(f"FATAL: Run {rid} missing from live ENA API!")
        all_reconciled = False
        continue

    e = ena_map[rid]
    md5_ena = e["fastq_md5"].split(";")
    ftp_ena = e["fastq_ftp"].split(";")

    md5_ok = (len(md5_ena) == 2) and (md5_ena[0] == s["fastq_1_md5"]) and (md5_ena[1] == s["fastq_2_md5"])
    ftp_ok = (len(ftp_ena) == 2) and (ftp_ena[0] == s["fastq_1_url"]) and (ftp_ena[1] == s["fastq_2_url"])
    alias_ok = e["sample_alias"] == s["sample_alias"]
    biosample_ok = e["sample_accession"] == s["biosample_accession"]
    experiment_ok = e["experiment_accession"] == s["experiment_accession"]
    reads_ok = int(e["read_count"]) == s["paired_reads_expected"]
    layout = e.get("library_layout", "PAIRED")

    if not (md5_ok and ftp_ok and alias_ok and biosample_ok and experiment_ok and reads_ok):
        all_reconciled = False
        print(f"DISCREPANCY in {rid}: md5={md5_ok}, ftp={ftp_ok}, alias={alias_ok}, biosample={biosample_ok}, exp={experiment_ok}, reads={reads_ok}")
    else:
        print(f"{rid:<13} {s['sample_alias']:<13} {s['group']:<10} {s['dose']:<4.0f} {s['dose_unit']:<3} {s['day']:<4} {layout:<8} {e['sample_accession']:<14} {e['experiment_accession']:<13} {int(e['read_count']):<12,d} MATCH_VERIFIED")

print("=" * 120)

if not all_reconciled:
    print("\nPILOT_MANIFEST_REVALIDATION_FAILED")
    raise SystemExit(1)
else:
    print("\nALL 12 PILOT RUNS 100% RECONCILED AGAINST LIVE ENA/SRA METADATA.")
