#!/usr/bin/env python3
"""Build and validate production human gene reference tables from authoritative Ensembl Release 112 and HGNC sources.

Architecture:
- Canonical gene universe originates from official, immutable Ensembl Release 112 GTF (Homo sapiens GRCh38.p14).
- Nomenclature, approved symbols, aliases, previous symbols, HGNC IDs, and Entrez IDs originate from a version-pinned HGNC snapshot.
- Ensembl genes lacking HGNC annotations are preserved in the canonical universe with unavailable nomenclature fields.
- Zero runtime hardcoding.
"""

import csv
import gzip
import hashlib
import json
import shutil
import sys
import time
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ENSEMBL_112_GTF_URL = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
ENSEMBL_112_EXPECTED_CHECKSUM = "8c87436bab973c871447743c2f5cbb2b557be0523fe6efe3ff33307a4c836002"
ENSEMBL_112_RELEASE_DATE = "2024-05"
GENOME_BUILD = "GRCh38.p14"
BUILDER_VERSION = "2.0.0"

HGNC_SNAPSHOT_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_SNAPSHOT_DATE = "2026-09-24"
HGNC_EXPECTED_CHECKSUM = "69bb5722d5a42bb355580deb2c9f197ce3d7f7d13807173b9674a7db65e52191"

HEADER = [
    "ensembl_gene_id",
    "gene_symbol",
    "entrez_gene_id",
    "hgnc_id",
    "gene_biotype",
    "synonyms",
    "chromosome",
    "description",
]


def sha256_file(path: Path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def ensure_file(url: str, dest_path: Path, expected_checksum: str | None = None) -> tuple[Path, str]:
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists():
        actual_cs = sha256_file(dest_path)
        if expected_checksum is None or actual_cs == expected_checksum:
            print(f"Using cached file at {dest_path} (SHA-256: {actual_cs[:12]}...)")
            return dest_path, actual_cs

    print(f"Downloading {url} to {dest_path}...")
    t0 = time.time()
    req = urllib.request.Request(url, headers={"User-Agent": "DNT-Organoid-AI/2.0"})
    with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as f_out:
        shutil.copyfileobj(resp, f_out)
    actual_cs = sha256_file(dest_path)
    print(f"Downloaded in {time.time() - t0:.2f}s ({dest_path.stat().st_size} bytes, SHA-256: {actual_cs})")
    if expected_checksum and actual_cs != expected_checksum:
        raise ValueError(f"Checksum mismatch for {dest_path}: expected {expected_checksum}, got {actual_cs}")
    return dest_path, actual_cs


def parse_ensembl_gtf(gtf_path: Path) -> dict[str, dict[str, str]]:
    """Parse canonical gene features from Ensembl 112 GTF."""
    print(f"Parsing canonical genes from {gtf_path}...")
    genes = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            
            chrom = parts[0]
            attrs = parts[8]
            info = {}
            for item in attrs.strip(";").split(";"):
                item = item.strip()
                if not item:
                    continue
                idx = item.find(" ")
                if idx != -1:
                    k = item[:idx]
                    v = item[idx+1:].strip('"')
                    info[k] = v
            gid = info.get("gene_id")
            if gid and gid.startswith("ENSG"):
                genes[gid] = {
                    "ensembl_gene_id": gid,
                    "ensembl_gene_name": info.get("gene_name", ""),
                    "gene_biotype": info.get("gene_biotype", "unknown"),
                    "chromosome": chrom,
                    "gene_version": info.get("gene_version", ""),
                }
    return genes


def parse_hgnc_snapshot(hgnc_path: Path) -> dict[str, dict[str, Any]]:
    """Parse nomenclature cross-reference layer from HGNC complete set snapshot."""
    print(f"Parsing nomenclature from {hgnc_path}...")
    hgnc_by_ensembl = {}
    opener = gzip.open if hgnc_path.name.endswith(".gz") else open
    with opener(hgnc_path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ens_id = (row.get("ensembl_gene_id") or "").strip()
            if not ens_id or not ens_id.startswith("ENSG"):
                continue
            
            symbol = (row.get("symbol") or "").strip()
            hgnc_id = (row.get("hgnc_id") or "").strip()
            name = (row.get("name") or "").strip()
            entrez = (row.get("entrez_id") or "").strip()
            
            raw_aliases = (row.get("alias_symbol") or "").strip()
            raw_prev = (row.get("prev_symbol") or "").strip()
            
            synonyms = set()
            for s in (raw_aliases.split("|") if raw_aliases else []) + (raw_prev.split("|") if raw_prev else []):
                s = s.strip().strip('"')
                if s and s != symbol:
                    synonyms.add(s)
                    no_hyphen = s.replace("-", "").replace(" ", "")
                    if no_hyphen and no_hyphen != symbol:
                        synonyms.add(no_hyphen)
            
            hgnc_by_ensembl[ens_id] = {
                "symbol": symbol,
                "hgnc_id": hgnc_id,
                "name": name,
                "entrez_id": entrez,
                "synonyms": sorted(synonyms),
            }
    return hgnc_by_ensembl


def build_production_reference(
    project_root: Path | None = None,
    output_path: Path | None = None,
) -> tuple[Path, str, dict[str, Any]]:
    """Build the production reference table cross-referencing Ensembl 112 with HGNC."""
    if project_root is None:
        project_root = Path(__file__).resolve().parents[1]
    
    raw_dir = project_root / "data" / "raw" / "references"
    gtf_path, gtf_checksum = ensure_file(
        ENSEMBL_112_GTF_URL,
        raw_dir / "Homo_sapiens.GRCh38.112.gtf.gz",
        ENSEMBL_112_EXPECTED_CHECKSUM,
    )

    # Check for version-controlled immutable archival HGNC snapshot
    archived_hgnc = project_root / "data" / "references" / "sources" / f"hgnc_complete_set_{HGNC_SNAPSHOT_DATE.replace('-', '')}.tsv.gz"
    if archived_hgnc.exists():
        hgnc_path = archived_hgnc
        uncomp = gzip.decompress(hgnc_path.read_bytes())
        hgnc_checksum = hashlib.sha256(uncomp).hexdigest()
        if hgnc_checksum != HGNC_EXPECTED_CHECKSUM:
            raise ValueError(f"Archived HGNC snapshot checksum mismatch: expected {HGNC_EXPECTED_CHECKSUM}, got {hgnc_checksum}")
        print(f"Using version-controlled archived HGNC snapshot at {hgnc_path} (uncompressed SHA-256: {hgnc_checksum[:12]}...)")
    else:
        hgnc_path, hgnc_checksum = ensure_file(
            HGNC_SNAPSHOT_URL,
            raw_dir / f"hgnc_complete_set_{HGNC_SNAPSHOT_DATE.replace('-', '')}.tsv",
            HGNC_EXPECTED_CHECKSUM,
        )
    
    ensembl_genes = parse_ensembl_gtf(gtf_path)
    hgnc_data = parse_hgnc_snapshot(hgnc_path)
    
    print("Combining Ensembl canonical genes with HGNC nomenclature...")
    records = []
    with_hgnc = 0
    with_symbol = 0
    with_hgnc_id = 0
    with_entrez = 0
    with_synonyms = 0
    biotype_counts = Counter()
    
    for gid, ens in ensembl_genes.items():
        h = hgnc_data.get(gid)
        biotype = ens["gene_biotype"]
        biotype_counts[biotype] += 1
        chrom = ens["chromosome"]
        
        if h:
            with_hgnc += 1
            sym = h["symbol"]
            h_id = h["hgnc_id"]
            ent = h["entrez_id"]
            desc = h["name"]
            syn_list = h["synonyms"]
            
            if sym:
                with_symbol += 1
            if h_id:
                with_hgnc_id += 1
            if ent:
                with_entrez += 1
            if syn_list:
                with_synonyms += 1
            
            syn_str = "|".join(syn_list)
        else:
            sym = ""
            h_id = ""
            ent = ""
            desc = ""
            syn_str = ""
        
        records.append((
            gid,
            sym,
            ent,
            h_id,
            biotype,
            syn_str,
            chrom,
            desc,
        ))
    
    # Deterministic sort by Ensembl ID
    records.sort(key=lambda r: r[0])
    
    # Write output TSV
    if output_path is None:
        output_path = project_root / "data" / "references" / "ensembl_human_genes_v112.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Writing {len(records)} records to {output_path}...")
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(records)
    
    output_checksum = sha256_file(output_path)
    print(f"Wrote reference table. SHA-256: {output_checksum}")
    
    # Write provenance JSON
    provenance = {
        "builder_version": BUILDER_VERSION,
        "generation_timestamp": datetime.now(timezone.utc).isoformat(),
        "canonical_gene_universe": {
            "source_authority": "Ensembl / EMBL-EBI",
            "release": "112",
            "genome_build": GENOME_BUILD,
            "release_date": ENSEMBL_112_RELEASE_DATE,
            "source_url": ENSEMBL_112_GTF_URL,
            "source_checksum_sha256": gtf_checksum,
            "total_canonical_genes": len(ensembl_genes),
        },
        "nomenclature_cross_reference": {
            "source_authority": "HUGO Gene Nomenclature Committee (HGNC)",
            "snapshot_date": HGNC_SNAPSHOT_DATE,
            "source_url": HGNC_SNAPSHOT_URL,
            "source_checksum_sha256": hgnc_checksum,
            "records_with_ensembl_ids": len(hgnc_data),
        },
        "coverage_metrics": {
            "total_canonical_genes": len(ensembl_genes),
            "genes_with_hgnc_cross_reference": with_hgnc,
            "genes_without_hgnc_cross_reference": len(ensembl_genes) - with_hgnc,
            "genes_with_approved_symbol": with_symbol,
            "genes_with_hgnc_id": with_hgnc_id,
            "genes_with_entrez_id": with_entrez,
            "genes_with_synonyms": with_synonyms,
            "top_biotypes": dict(biotype_counts.most_common(10)),
        },
        "production_reference_asset": {
            "file_path": str(output_path.relative_to(project_root)).replace("\\", "/"),
            "file_size_bytes": output_path.stat().st_size,
            "record_count": len(records),
            "sha256_checksum": output_checksum,
        },
    }
    
    prov_path = output_path.parent / "ensembl_human_genes_v112_provenance.json"
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2)
    print(f"Wrote provenance metadata to {prov_path}")
    
    return output_path, output_checksum, provenance


if __name__ == "__main__":
    out_file, cs, prov = build_production_reference()
    print("\n=== Reference Build Summary ===")
    print(f"Total Canonical Ensembl 112 Genes: {prov['coverage_metrics']['total_canonical_genes']}")
    print(f"Genes with HGNC symbols: {prov['coverage_metrics']['genes_with_approved_symbol']}")
    print(f"Genes without HGNC symbols: {prov['coverage_metrics']['genes_without_hgnc_cross_reference']}")
    print(f"Production TSV Checksum: {cs}")
