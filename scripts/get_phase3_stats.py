"""Print canonical metrics for Phase 3 Scientific Verification."""

import json
import pathlib
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
METADATA_DIR = REPO_ROOT / "data" / "metadata"
PHASE3_DIR = REPO_ROOT / "outputs" / "verification" / "phase3"

df_inv = pd.read_csv(METADATA_DIR / "dnt_master_dataset_inventory.csv")
df_samples = pd.read_csv(METADATA_DIR / "dnt_sample_metadata_v1.csv")
df_cond = pd.read_csv(METADATA_DIR / "dnt_experimental_conditions_inventory.csv")

print("=" * 80)
print("PHASE 3 VERIFIED SCIENTIFIC METRICS")
print("=" * 80)

print("\n--- 1. DATASET ELIGIBILITY (1,666 CENSUS) ---")
for dec, cnt in df_inv["decision"].value_counts().items():
    print(f"  {dec:<20}: {cnt:>5} ({cnt/len(df_inv)*100:>5.2f}%)")
print(f"  TOTAL               : {len(df_inv):>5} (100.00%)")

included_inv = df_inv[df_inv["decision"] == "INCLUDE"]
organoids = [m for m in included_inv["biological_model"] if any(k in m.lower() for k in ["organoid", "spheroid", "assembloid"])]
other_neural = [m for m in included_inv["biological_model"] if m not in organoids]

print("\n--- 2. INCLUDED BIOLOGICAL MODELS (48 DATASETS) ---")
print(f"  Human Brain Organoids / Spheroids / Assembloids : {len(organoids):>3} datasets")
print(f"  Other Neural Models (NPCs, NSCs, Neurons, LUHMES): {len(other_neural):>3} datasets")
for model, count in included_inv["biological_model"].value_counts().items():
    print(f"    - {model:<45}: {count:>2} datasets")

print("\n--- 3. SAMPLE-LEVEL METADATA (1,298 RNA-seq SAMPLES) ---")
print(f"  Total Verified RNA-seq Samples: {len(df_samples):>5}")
n_treat = sum(df_samples["Treatment / Control Group"] == "Treatment")
n_ctrl = sum(df_samples["Treatment / Control Group"] == "Control")
print(f"  Chemical Treatment Samples    : {n_treat:>5} ({n_treat/len(df_samples)*100:.2f}%)")
print(f"  Matched Vehicle Controls      : {n_ctrl:>5} ({n_ctrl/len(df_samples)*100:.2f}%)")
for lbl, count in df_samples["DNT Outcome Label"].value_counts().items():
    print(f"    - Outcome Label [{lbl}]: {count}")

print("\n--- 4. EXPERIMENTAL CONDITIONS INVENTORY (219 CONDITIONS) ---")
print(f"  Total Distinct Conditions      : {len(df_cond):>5}")
n_adeq = sum(df_cond["Replication Status"].str.startswith("ADEQUATE"))
n_single = sum(df_cond["Replication Status"].str.startswith("SINGLE"))
print(f"  Adequately Replicated (n >= 2) : {n_adeq:>5} ({n_adeq/len(df_cond)*100:.2f}%)")
print(f"  Single Replicate (Flagged)     : {n_single:>5} ({n_single/len(df_cond)*100:.2f}%)")
n_matched_ctrl = sum(df_cond["Control Status"] == "MATCHED_VEHICLE_CONTROL")
print(f"  Matched Vehicle Controls       : {n_matched_ctrl:>5} ({n_matched_ctrl/len(df_cond)*100:.2f}%)")
print(f"  Unique Verified Chemicals      : {df_cond['Compound Name'].nunique():>5}")
print(f"  Chemical Partition Groups      : {df_cond['Chemical Partition Group'].nunique():>5}")

print("\n--- 5. METADATA COMPLETENESS (1,298 SAMPLES) ---")
with open(PHASE3_DIR / "phase3_metadata_completeness.json", "r", encoding="utf-8") as f:
    comp = json.load(f)
for field in [
    "Dataset ID", "Sample ID", "Species", "Cell Type / Organoid Type", "Cell Line / Donor",
    "Organoid Age / Differentiation Day", "Age at Exposure", "Age at RNA Collection",
    "Compound Name", "Compound ID", "Compound Concentration", "Exposure Duration",
    "Exposure Pattern", "Treatment / Control Group", "Biological Replicates",
    "RNA-seq Type", "Sequencing Platform", "Study / Batch ID", "Gene Expression Data",
    "DNT Outcome Label", "DNT Evidence Source", "DNT Label Confidence"
]:
    data = comp.get(field, {})
    vc = data.get("valid_count", 0)
    pct = data.get("actual_completeness_pct", 0.0)
    print(f"  {field:<36}: {vc:>4}/1298 ({pct:>6.2f}%)")
