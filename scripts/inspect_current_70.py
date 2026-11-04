#!/usr/bin/env python3
"""Inspect and audit the 70 datasets currently classified as INCLUDE in Phase 1."""

import json
import pathlib
import sys
import pandas as pd

if sys.stdout.encoding != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
INVENTORY_CSV = REPO_ROOT / "data" / "metadata" / "dnt_master_dataset_inventory.csv"
CACHE_DIR = REPO_ROOT / "outputs" / "verification" / "cache"
GDS_CACHE_FILE = CACHE_DIR / "gds_esummary_cache.json"
PUBMED_CACHE_FILE = CACHE_DIR / "pubmed_cache.json"

df = pd.read_csv(INVENTORY_CSV)
inc = df[df["decision"] == "INCLUDE"]

with open(GDS_CACHE_FILE, "r", encoding="utf-8") as f:
    gds_cache = json.load(f)

print(f"Total datasets classified as INCLUDE: {len(inc)}")
print("=" * 100)

for i, (_, r) in enumerate(inc.iterrows(), 1):
    gse = r["gse"]
    gds = gds_cache.get(gse, {})
    title = gds.get("title", r.get("title", ""))
    summary = gds.get("summary", "")
    n_samples = r["n_samples"]
    model = r["biological_model"]
    compounds = r["compounds"]
    pmids = r["pubmed_ids"]
    
    print(f"[{i:02d}/70] {gse} (n={n_samples}) | Model: {model} | Reported Comps: {compounds}")
    print(f"     Title: {title[:90]}")
    # Show first 200 chars of summary
    print(f"     Summary excerpt: {summary[:180]}...")
    print("-" * 100)
