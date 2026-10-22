"""Verify authoritative ENA/SRA metadata for the 12 approved pilot runs."""

import csv
import urllib.request

from src.quantification_bridge.manifest import APPROVED_PILOT_RUNS

print("Re-querying ENA API for PRJNA802072 runs...")
url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA802072&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,sample_alias,sample_title,fastq_ftp,fastq_md5,read_count,base_count&format=tsv&download=true"

req = urllib.request.Request(url, headers={"User-Agent": "DNT-Organoid-AI/2.0"})
with urllib.request.urlopen(req) as resp:
    lines = resp.read().decode("utf-8").strip().split("\n")

reader = csv.DictReader(lines, delimiter="\t")
ena_runs = {r["run_accession"]: r for r in reader}
print(f"Total runs returned by ENA API for PRJNA802072: {len(ena_runs)}")

# Verify each of the 12 approved pilot runs
discrepancies = []
for pilot in APPROVED_PILOT_RUNS:
    run_id = pilot.run_accession
    if run_id not in ena_runs:
        discrepancies.append(f"Run {run_id} missing from ENA API response!")
        continue

    rec = ena_runs[run_id]
    if rec["sample_accession"] != pilot.biosample_accession:
        discrepancies.append(f"Run {run_id} BioSample mismatch: expected {pilot.biosample_accession}, got {rec['sample_accession']}")
    if rec["experiment_accession"] != pilot.experiment_accession:
        discrepancies.append(f"Run {run_id} Experiment mismatch: expected {pilot.experiment_accession}, got {rec['experiment_accession']}")
    # In ENA, sample_alias holds the original sample title (e.g. VPA_0_1_1), while sample_title has generic text
    if rec["sample_alias"] != pilot.sample_title:
        discrepancies.append(f"Run {run_id} Sample alias mismatch: expected {pilot.sample_title}, got {rec['sample_alias']}")
    if int(rec["read_count"]) != pilot.paired_reads_expected:
        discrepancies.append(f"Run {run_id} Read count mismatch: expected {pilot.paired_reads_expected}, got {rec['read_count']}")
    if int(rec["base_count"]) != pilot.total_bases_expected:
        discrepancies.append(f"Run {run_id} Base count mismatch: expected {pilot.total_bases_expected}, got {rec['base_count']}")

    # Check FTP and MD5s
    ftp_parts = rec["fastq_ftp"].split(";")
    md5_parts = rec["fastq_md5"].split(";")
    if len(ftp_parts) != 2 or len(md5_parts) != 2:
        discrepancies.append(f"Run {run_id} does not have exactly 2 FASTQ files in ENA!")
    else:
        if ftp_parts[0] != pilot.fastq_1_url or ftp_parts[1] != pilot.fastq_2_url:
            discrepancies.append(f"Run {run_id} FTP URL mismatch")
        if md5_parts[0] != pilot.fastq_1_md5 or md5_parts[1] != pilot.fastq_2_md5:
            discrepancies.append(f"Run {run_id} MD5 mismatch: expected {pilot.fastq_1_md5};{pilot.fastq_2_md5}, got {rec['fastq_md5']}")

if discrepancies:
    print("METADATA VERIFICATION FAILED:")
    for d in discrepancies:
        print("  -", d)
    raise SystemExit(1)
else:
    print("ALL 12 PILOT RUNS 100% VERIFIED AGAINST LIVE ENA/SRA API!")
    print("BioProject: PRJNA802072 | Study: SRP357550 | Contrast: Day 28 VPA 250 uM vs 0 uM (N=6 vs N=6)")
    print("All BioSample accessions, Experiment accessions, Sample Aliases, Read Counts, Base Counts, and MD5s match.")
