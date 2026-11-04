import json
import pandas as pd

inv = pd.read_csv('data/metadata/dnt_master_dataset_inventory.csv')
samples = pd.read_csv('data/metadata/dnt_sample_metadata_v1.csv')
cond = pd.read_csv('data/metadata/dnt_experimental_conditions_inventory.csv')
with open('outputs/verification/phase3/phase3_metadata_completeness.json') as f:
    comp = json.load(f)

print("=== DATASET DECISIONS (TOTAL 1666) ===")
print(inv['decision'].value_counts())

print("\n=== BIOLOGICAL MODELS (INCLUDED) ===")
inc = inv[inv['decision'] == 'INCLUDE']
print(inc['biological_model'].value_counts())

print("\n=== BIOLOGICAL MODELS (REFERENCE DATASETS) ===")
ref = inv[inv['decision'] == 'REFERENCE_DATASET']
print(ref['biological_model'].value_counts())

print("\n=== SAMPLES SUMMARY ===")
print(f"Total samples: {len(samples)}")
print(samples['Treatment / Control Group'].value_counts())
treat_samples = samples[samples['Treatment / Control Group'] == 'Treatment']
unique_chems = treat_samples['Compound Name'].nunique()
print(f"Unique verified chemicals: {unique_chems}")
print(f"Unique datasets with samples: {samples['Dataset ID'].nunique()}")

print("\n=== EXPERIMENTAL CONDITIONS SUMMARY ===")
print(f"Total experimental conditions: {len(cond)}")
print("Replication status:")
print(cond['Replication Status'].value_counts())
print("\nControl status:")
print(cond['Control Status'].value_counts())
print(f"Unique chemicals in conditions: {cond['Compound Name'].nunique()}")
print(f"Unique partition groups: {cond['Chemical Partition Group'].nunique()}")
print(f"Unique studies in conditions: {cond['Dataset ID'].nunique()}")

print("\n=== METADATA COMPLETENESS ===")
for k in ["Compound Concentration", "Exposure Duration", "Organoid Age / Differentiation Day", "Age at Exposure", "Age at RNA Collection"]:
    v = comp.get(k, {})
    print(f"{k:35}: {v.get('valid_count')} / {v.get('total_samples')} ({v.get('actual_completeness_pct')}%)")

print("\n=== RECOVERED FROM PREVIOUS EXCLUSIONS ===")
recovered = inv[(inv['phase1_decision'] == 'EXCLUDE') & (inv['decision'] == 'INCLUDE')]
print(f"Recovered from Phase 1 EXCLUDE to INCLUDE: {len(recovered)}")
print(recovered[['gse', 'biological_model', 'compounds']].to_string())

print("\n=== RECOVERED FROM NEEDS_REVIEW ===")
rec_nr = inv[(inv['phase1_decision'] == 'NEEDS_REVIEW') & (inv['decision'] == 'INCLUDE')]
print(f"Recovered from Phase 1 NEEDS_REVIEW to INCLUDE: {len(rec_nr)}")
print(rec_nr[['gse', 'biological_model', 'compounds']].to_string())
