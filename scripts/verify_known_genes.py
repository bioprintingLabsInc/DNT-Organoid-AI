import sys
from pathlib import Path

# Add project root to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.gene_harmonization import (
    GeneIdentifierType,
    GeneReference,
    load_config,
)

def verify():
    config = load_config()
    ref = GeneReference.load_from_config(config)

    queries = [
        ("SOX2", GeneIdentifierType.GENE_SYMBOL),
        ("PAX6", GeneIdentifierType.GENE_SYMBOL),
        ("FOXG1", GeneIdentifierType.GENE_SYMBOL),
        ("MAP2", GeneIdentifierType.GENE_SYMBOL),
        ("SYN1", GeneIdentifierType.GENE_SYMBOL),
        ("HMOX1", GeneIdentifierType.GENE_SYMBOL),
        ("CASP3", GeneIdentifierType.GENE_SYMBOL),
        ("MKI67", GeneIdentifierType.GENE_SYMBOL),
        ("KI67", GeneIdentifierType.GENE_SYMBOL),
        ("Ki-67", GeneIdentifierType.GENE_SYMBOL),
        ("EOMES", GeneIdentifierType.GENE_SYMBOL),
        ("TBR2", GeneIdentifierType.GENE_SYMBOL),
        ("ENSG00000181449", GeneIdentifierType.ENSEMBL_GENE_ID),   # SOX2
        ("ENSG00000007372.13", GeneIdentifierType.ENSEMBL_GENE_ID), # PAX6 versioned
        ("ENSG00000228037", GeneIdentifierType.ENSEMBL_GENE_ID),   # Ensembl 112 gene WITHOUT HGNC annotation
        ("FLIP", GeneIdentifierType.GENE_SYMBOL),                  # Ambiguous synonym
    ]

    header = f"{'Query':<20} | {'Type':<16} | {'Canonical Ensembl ID':<19} | {'Approved Symbol':<16} | {'Status':<16} | {'Reference Provenance'}"
    print(header)
    print("-" * len(header))
    for q, t in queries:
        res = ref.lookup(q, t)
        provenance = f"{ref.reference_id}:v{ref.reference_version} ({ref.checksum[:8]}...)"
        print(f"{q:<20} | {t.value:<16} | {str(res.canonical_gene_id):<19} | {str(res.approved_symbol):<16} | {res.status.value:<16} | {provenance}")

if __name__ == "__main__":
    verify()
