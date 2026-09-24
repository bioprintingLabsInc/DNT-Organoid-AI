"""Versioned human gene reference catalog, loading, indexing, and lookup."""

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .models import GeneIdentifierType, GeneMappingStatus


@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    ensembl_gene_id: str
    gene_symbol: str
    entrez_gene_id: str
    hgnc_id: str
    gene_biotype: str
    synonyms: tuple[str, ...]
    chromosome: str
    description: str


@dataclass(frozen=True, slots=True)
class ReferenceLookupResult:
    status: GeneMappingStatus
    canonical_gene_id: str | None
    approved_symbol: str | None
    candidate_canonical_ids: tuple[str, ...]
    reason: str | None


class GeneReference:
    """Indexed human gene reference for deterministic identifier resolution."""

    def __init__(
        self,
        reference_id: str,
        reference_version: str,
        genome_build: str,
        source_authority: str,
        records: tuple[ReferenceRecord, ...],
        checksum: str,
    ) -> None:
        self.reference_id = reference_id
        self.reference_version = reference_version
        self.genome_build = genome_build
        self.source_authority = source_authority
        self.records = records
        self.checksum = checksum

        self._by_ensembl: dict[str, ReferenceRecord] = {}
        self._by_symbol: dict[str, list[ReferenceRecord]] = {}
        self._by_entrez: dict[str, list[ReferenceRecord]] = {}
        self._by_hgnc: dict[str, list[ReferenceRecord]] = {}
        self._by_synonym: dict[str, list[ReferenceRecord]] = {}

        self._build_indexes()

    def _build_indexes(self) -> None:
        for r in self.records:
            self._by_ensembl[r.ensembl_gene_id] = r
            if r.gene_symbol:
                self._by_symbol.setdefault(r.gene_symbol.upper(), []).append(r)
            if r.entrez_gene_id:
                self._by_entrez.setdefault(r.entrez_gene_id, []).append(r)
            if r.hgnc_id:
                self._by_hgnc.setdefault(r.hgnc_id.upper(), []).append(r)
            for syn in r.synonyms:
                if syn:
                    self._by_synonym.setdefault(syn.upper(), []).append(r)

    @classmethod
    def load_from_config(cls, config: dict[str, Any], project_root: Path | None = None) -> "GeneReference":
        """Load reference table using the specification in the configuration."""
        ref_spec = config.get("reference", {})
        rel_path = ref_spec.get("file_path", "data/references/ensembl_human_genes_v112.tsv")
        if project_root is None:
            project_root = Path(__file__).resolve().parents[2]
        file_path = project_root / rel_path if not Path(rel_path).is_absolute() else Path(rel_path)

        if not file_path.exists():
            raise FileNotFoundError(f"Gene reference file not found at: {file_path}")

        content_bytes = file_path.read_bytes()
        calculated_checksum = hashlib.sha256(content_bytes).hexdigest()

        text = content_bytes.decode("utf-8-sig")
        reader = csv.reader(text.splitlines(), delimiter="\t")
        rows = list(reader)
        if not rows:
            raise ValueError(f"Gene reference file is empty: {file_path}")

        header = [col.strip() for col in rows[0]]
        expected_cols = {"ensembl_gene_id", "gene_symbol"}
        if not expected_cols.issubset(set(header)):
            raise ValueError(f"Gene reference header missing required columns {expected_cols}: {header}")

        col_idx = {name: i for i, name in enumerate(header)}
        records = []
        for line_num, row in enumerate(rows[1:], start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"Row {line_num} in reference has width {len(row)}, expected {len(header)}")

            def get_val(col_name: str) -> str:
                idx = col_idx.get(col_name)
                return row[idx].strip() if idx is not None and idx < len(row) else ""

            ens_id = get_val("ensembl_gene_id")
            symbol = get_val("gene_symbol")
            entrez = get_val("entrez_gene_id")
            hgnc = get_val("hgnc_id")
            biotype = get_val("gene_biotype")
            synonyms_str = get_val("synonyms")
            synonyms = tuple(s.strip() for s in synonyms_str.split("|") if s.strip()) if synonyms_str else ()
            chrom = get_val("chromosome")
            desc = get_val("description")

            records.append(
                ReferenceRecord(
                    ensembl_gene_id=ens_id,
                    gene_symbol=symbol,
                    entrez_gene_id=entrez,
                    hgnc_id=hgnc,
                    gene_biotype=biotype,
                    synonyms=synonyms,
                    chromosome=chrom,
                    description=desc,
                )
            )

        return cls(
            reference_id=ref_spec.get("reference_id", "ensembl_human_genes"),
            reference_version=ref_spec.get("reference_version", "112"),
            genome_build=ref_spec.get("genome_build", "GRCh38.p14"),
            source_authority=ref_spec.get("source_authority", "Ensembl / EMBL-EBI"),
            records=tuple(records),
            checksum=calculated_checksum,
        )

    def lookup(self, raw_id: str, id_type: GeneIdentifierType) -> ReferenceLookupResult:
        """Query reference using the declared/detected identifier type."""
        stripped = raw_id.strip()
        if not stripped:
            return ReferenceLookupResult(
                status=GeneMappingStatus.INVALID,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=(),
                reason="Identifier is empty or whitespace.",
            )

        if id_type == GeneIdentifierType.ENSEMBL_GENE_ID:
            return self.lookup_ensembl(stripped)
        elif id_type == GeneIdentifierType.GENE_SYMBOL:
            return self.lookup_symbol(stripped)
        elif id_type == GeneIdentifierType.ENTREZ_GENE_ID:
            return self.lookup_entrez(stripped)
        elif id_type == GeneIdentifierType.HGNC_ID:
            return self.lookup_hgnc(stripped)
        else:
            return ReferenceLookupResult(
                status=GeneMappingStatus.INVALID,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=(),
                reason=f"Unsupported or unconfigured identifier type: '{id_type}'.",
            )

    def lookup_ensembl(self, raw_id: str) -> ReferenceLookupResult:
        """Resolve Ensembl Gene ID, stripping version suffix while keeping original unchanged."""
        base_id = raw_id.split(".")[0].strip()
        record = self._by_ensembl.get(base_id)
        if record is not None:
            return ReferenceLookupResult(
                status=GeneMappingStatus.UNIQUELY_MAPPED,
                canonical_gene_id=record.ensembl_gene_id,
                approved_symbol=record.gene_symbol or None,
                candidate_canonical_ids=(record.ensembl_gene_id,),
                reason=None,
            )
        return ReferenceLookupResult(
            status=GeneMappingStatus.UNMAPPED,
            canonical_gene_id=None,
            approved_symbol=None,
            candidate_canonical_ids=(),
            reason=f"Ensembl Gene ID '{raw_id}' not found in reference '{self.reference_id}' version '{self.reference_version}'.",
        )

    def lookup_symbol(self, raw_symbol: str) -> ReferenceLookupResult:
        """Resolve gene symbol using approved symbols, then synonyms, preserving ambiguity."""
        key = raw_symbol.upper()
        hits = self._by_symbol.get(key, [])
        if len(hits) == 1:
            rec = hits[0]
            return ReferenceLookupResult(
                status=GeneMappingStatus.UNIQUELY_MAPPED,
                canonical_gene_id=rec.ensembl_gene_id,
                approved_symbol=rec.gene_symbol,
                candidate_canonical_ids=(rec.ensembl_gene_id,),
                reason=None,
            )
        elif len(hits) > 1:
            candidate_ids = tuple(sorted({h.ensembl_gene_id for h in hits}))
            return ReferenceLookupResult(
                status=GeneMappingStatus.AMBIGUOUS,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=candidate_ids,
                reason=f"Gene symbol '{raw_symbol}' matches multiple canonical Ensembl IDs: {', '.join(candidate_ids)}.",
            )

        # Check synonyms
        syn_hits = self._by_synonym.get(key, [])
        if len(syn_hits) == 1:
            rec = syn_hits[0]
            return ReferenceLookupResult(
                status=GeneMappingStatus.UNIQUELY_MAPPED,
                canonical_gene_id=rec.ensembl_gene_id,
                approved_symbol=rec.gene_symbol,
                candidate_canonical_ids=(rec.ensembl_gene_id,),
                reason=None,
            )
        elif len(syn_hits) > 1:
            candidate_ids = tuple(sorted({h.ensembl_gene_id for h in syn_hits}))
            return ReferenceLookupResult(
                status=GeneMappingStatus.AMBIGUOUS,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=candidate_ids,
                reason=f"Gene synonym '{raw_symbol}' matches multiple canonical Ensembl IDs: {', '.join(candidate_ids)}.",
            )

        return ReferenceLookupResult(
            status=GeneMappingStatus.UNMAPPED,
            canonical_gene_id=None,
            approved_symbol=None,
            candidate_canonical_ids=(),
            reason=f"Gene symbol '{raw_symbol}' not found in reference '{self.reference_id}' version '{self.reference_version}'.",
        )

    def lookup_entrez(self, raw_entrez: str) -> ReferenceLookupResult:
        """Resolve Entrez / NCBI Gene ID."""
        hits = self._by_entrez.get(raw_entrez, [])
        if len(hits) == 1:
            rec = hits[0]
            return ReferenceLookupResult(
                status=GeneMappingStatus.UNIQUELY_MAPPED,
                canonical_gene_id=rec.ensembl_gene_id,
                approved_symbol=rec.gene_symbol or None,
                candidate_canonical_ids=(rec.ensembl_gene_id,),
                reason=None,
            )
        elif len(hits) > 1:
            candidate_ids = tuple(sorted({h.ensembl_gene_id for h in hits}))
            return ReferenceLookupResult(
                status=GeneMappingStatus.AMBIGUOUS,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=candidate_ids,
                reason=f"Entrez Gene ID '{raw_entrez}' matches multiple canonical Ensembl IDs: {', '.join(candidate_ids)}.",
            )
        return ReferenceLookupResult(
            status=GeneMappingStatus.UNMAPPED,
            canonical_gene_id=None,
            approved_symbol=None,
            candidate_canonical_ids=(),
            reason=f"Entrez Gene ID '{raw_entrez}' not found in reference '{self.reference_id}' version '{self.reference_version}'.",
        )

    def lookup_hgnc(self, raw_hgnc: str) -> ReferenceLookupResult:
        """Resolve HGNC ID."""
        key = raw_hgnc.upper()
        hits = self._by_hgnc.get(key, [])
        if len(hits) == 1:
            rec = hits[0]
            return ReferenceLookupResult(
                status=GeneMappingStatus.UNIQUELY_MAPPED,
                canonical_gene_id=rec.ensembl_gene_id,
                approved_symbol=rec.gene_symbol or None,
                candidate_canonical_ids=(rec.ensembl_gene_id,),
                reason=None,
            )
        elif len(hits) > 1:
            candidate_ids = tuple(sorted({h.ensembl_gene_id for h in hits}))
            return ReferenceLookupResult(
                status=GeneMappingStatus.AMBIGUOUS,
                canonical_gene_id=None,
                approved_symbol=None,
                candidate_canonical_ids=candidate_ids,
                reason=f"HGNC ID '{raw_hgnc}' matches multiple canonical Ensembl IDs: {', '.join(candidate_ids)}.",
            )
        return ReferenceLookupResult(
            status=GeneMappingStatus.UNMAPPED,
            canonical_gene_id=None,
            approved_symbol=None,
            candidate_canonical_ids=(),
            reason=f"HGNC ID '{raw_hgnc}' not found in reference '{self.reference_id}' version '{self.reference_version}'.",
        )
