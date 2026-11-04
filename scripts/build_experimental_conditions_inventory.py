"""Generate canonical DNT Experimental Conditions Inventory.

Preserves every distinct experimental condition across all included datasets:
- Groups by (Dataset ID, Compound Name, Concentration, Duration, Pattern, Developmental Stage, Model, Cell Line)
- Links sample GSM accessions, replicate counts, and matched vehicle control GSMs
- Normalizes concentrations to micromolar (uM) and durations to hours where scientifically valid
- Assigns chemical partition groups (CHEM_{slug}) ensuring cross-validation isolation
- Flags single-replicate or unmatched conditions for review
- Generates CSV and Excel formats in data/metadata/ and outputs/verification/phase3/
"""

import json
import pathlib
import re
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "data" / "metadata"
PHASE3_DIR = REPO_ROOT / "outputs" / "verification" / "phase3"

def normalize_conc(dose_str: str) -> str:
    s = str(dose_str).strip()
    if s in ["0", "0.0", "control", "none", "untreated"]:
        return "0.0"
    m = re.match(r"^([\d\.]+)\s*([a-zA-Zµ/]+)$", s)
    if not m:
        return "NOT_MAPPED"
    val = float(m.group(1))
    unit = m.group(2).lower()
    if unit in ["mm", "millimolar"]:
        return str(round(val * 1000.0, 4))
    elif unit in ["um", "µm", "micromolar"]:
        return str(round(val, 4))
    elif unit in ["nm", "nanomolar"]:
        return str(round(val * 0.001, 6))
    return "NOT_MAPPED"

def normalize_dur(dur_str: str) -> str:
    s = str(dur_str).strip()
    m = re.match(r"^([\d\.]+)\s*([a-zA-Z]+)$", s)
    if not m:
        return "NOT_REPORTED"
    val = float(m.group(1))
    unit = m.group(2).lower()
    if unit in ["hours", "hour", "hrs", "hr", "h"]:
        return str(round(val, 2))
    elif unit in ["days", "day", "d"]:
        return str(round(val * 24.0, 2))
    elif unit in ["weeks", "week", "w"]:
        return str(round(val * 168.0, 2))
    elif unit in ["minutes", "minute", "min"]:
        return str(round(val / 60.0, 2))
    return "NOT_REPORTED"

def main():
    print("=" * 100)
    print("BUILDING DNT EXPERIMENTAL CONDITIONS INVENTORY")
    print("=" * 100)

    df_samples = pd.read_csv(METADATA_DIR / "dnt_sample_metadata_v1.csv")
    print(f"Loaded {len(df_samples)} verified sample records across {df_samples['Dataset ID'].nunique()} datasets.")

    # 1. Map controls by GSE
    gse_controls = {}
    for _, s in df_samples.iterrows():
        if s["Treatment / Control Group"] == "Control":
            gse_controls.setdefault(s["Dataset ID"], []).append(s["Sample ID"])

    # 2. Cluster treatment samples into distinct conditions
    cond_map = {}
    for _, s in df_samples.iterrows():
        if s["Treatment / Control Group"] == "Treatment":
            key = (
                str(s["Dataset ID"]).strip(),
                str(s["Study / Batch ID"]).strip(),
                str(s["Compound Name"]).strip(),
                str(s["Compound Concentration"]).strip(),
                str(s["Exposure Duration"]).strip(),
                str(s["Exposure Pattern"]).strip(),
                str(s["Organoid Age / Differentiation Day"]).strip(),
                str(s["Cell Type / Organoid Type"]).strip(),
                str(s["Cell Line / Donor"]).strip(),
            )
            cond_map.setdefault(key, []).append(s["Sample ID"])

    print(f"Identified {len(cond_map)} distinct experimental conditions.")

    rows = []
    gse_cond_counters = {}
    for key, sample_ids in cond_map.items():
        gse, study_id, chem_name, dose, dur, pattern, stage, model, cell_line = key
        gse_cond_counters[gse] = gse_cond_counters.get(gse, 0) + 1
        cond_idx = gse_cond_counters[gse]

        n_reps = len(sample_ids)
        ctrl_gsms = gse_controls.get(gse, [])

        rep_status = "ADEQUATE_REPLICATION (n>=2)" if n_reps >= 2 else "SINGLE_REPLICATE (FLAGGED_FOR_REVIEW)"
        ctrl_status = "MATCHED_VEHICLE_CONTROL" if len(ctrl_gsms) >= 1 else "UNMATCHED_FLAGGED_FOR_REVIEW"
        eligibility = "INCLUDE" if (n_reps >= 2 and len(ctrl_gsms) >= 1) else "NEEDS_REVIEW_CONDITION"

        # Safe partition grouping slug
        clean_chem = re.sub(r"[^a-zA-Z0-9]+", "_", chem_name.strip()).strip("_").upper()
        partition_group = f"CHEM_{clean_chem}"

        norm_c = normalize_conc(dose)
        norm_d = normalize_dur(dur)

        evidence = f"NCBI GEO {gse} ({len(sample_ids)} treated samples: {sample_ids[0]}..{sample_ids[-1]}; matched controls: {len(ctrl_gsms)} samples)"

        rows.append({
            "Condition ID": f"COND_{gse}_{cond_idx:03d}",
            "Dataset ID": gse,
            "Study / BioProject ID": study_id,
            "Compound Name": chem_name,
            "Chemical Partition Group": partition_group,
            "Concentration Raw": dose,
            "Concentration Normalized (uM)": norm_c,
            "Exposure Duration Raw": dur,
            "Exposure Duration (Hours)": norm_d,
            "Exposure Pattern": pattern,
            "Developmental Stage": stage,
            "Biological Model": model,
            "Cell Line / Donor": cell_line,
            "Sample GSMs": ", ".join(sample_ids),
            "Replicate Count (n)": n_reps,
            "Replication Status": rep_status,
            "Matched Control GSMs": ", ".join(ctrl_gsms) if ctrl_gsms else "NOT_REPORTED",
            "Control Status": ctrl_status,
            "Condition Eligibility": eligibility,
            "Evidence Source": evidence,
        })

    df_conditions = pd.DataFrame(rows)
    
    # Save outputs
    csv_meta = METADATA_DIR / "dnt_experimental_conditions_inventory.csv"
    xlsx_meta = METADATA_DIR / "dnt_experimental_conditions_inventory.xlsx"
    csv_p3 = PHASE3_DIR / "dnt_experimental_conditions_inventory.csv"
    xlsx_p3 = PHASE3_DIR / "dnt_experimental_conditions_inventory.xlsx"

    df_conditions.to_csv(csv_meta, index=False)
    df_conditions.to_excel(xlsx_meta, index=False)
    df_conditions.to_csv(csv_p3, index=False)
    df_conditions.to_excel(xlsx_p3, index=False)

    print(f"Experimental conditions inventory successfully saved to:\n  {csv_meta}\n  {xlsx_meta}\n  {csv_p3}\n  {xlsx_p3}")
    print(f"Total Conditions: {len(df_conditions)}")
    print(f"Adequately Replicated Conditions (n>=2): {sum(df_conditions['Replicate Count (n)'] >= 2)} ({sum(df_conditions['Replicate Count (n)'] >= 2)/len(df_conditions)*100:.1f}%)")
    print(f"Single Replicate Conditions (Flagged for Review): {sum(df_conditions['Replicate Count (n)'] == 1)}")
    print(f"Matched Control Rate: {sum(df_conditions['Control Status'] == 'MATCHED_VEHICLE_CONTROL') / len(df_conditions)*100:.1f}%")
    print(f"Unique Partition Groups: {df_conditions['Chemical Partition Group'].nunique()}")

if __name__ == "__main__":
    main()
