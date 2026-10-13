"""Governed gene mapping adapter resolving Ensembl canonical IDs to Entrez Gene IDs for Stage 8C."""

import math
from pathlib import Path
from typing import Sequence

from src.differential_expression.models import GeneDifferentialExpressionResult
from src.gene_harmonization.reference import GeneReference

from .models import (
    MappingStatus,
    PathwayGeneMappingRecord,
    PathwayGeneMappingResult,
)


class EnsemblToEntrezAdapter:
    """Consumer-side mapping adapter resolving canonical Ensembl 112 genes to NCBI Entrez IDs."""

    def __init__(self, gene_reference: GeneReference | None = None) -> None:
        if gene_reference is None:
            config = {"reference": {"file_path": "data/references/ensembl_human_genes_v112.tsv"}}
            project_root = Path(__file__).resolve().parents[3]
            self.gene_reference = GeneReference.load_from_config(config, project_root=project_root)
        else:
            self.gene_reference = gene_reference

    def map_genes(
        self, gene_results: Sequence[GeneDifferentialExpressionResult]
    ) -> PathwayGeneMappingResult:
        """Map evaluated Stage 7 gene results to Entrez Gene IDs using governed local reference."""
        records: list[PathwayGeneMappingRecord] = []
        entrez_counts: dict[str, int] = {}
        mapped_count = 0
        unmapped_count = 0
        invalid_stat_count = 0
        one_to_many_count = 0

        # First pass: identify mapping and count Entrez occurrences to detect duplicates
        intermediate: list[tuple[GeneDifferentialExpressionResult, str | None, str, str | None]] = []

        for gene in gene_results:
            canonical_id = gene.canonical_gene_id
            ref_record = self.gene_reference._by_ensembl.get(canonical_id)

            # Check statistical validity
            is_stat_valid = (
                gene.result_status == "OK"
                and gene.wald_statistic is not None
                and not math.isnan(gene.wald_statistic)
                and not math.isinf(gene.wald_statistic)
            )

            if ref_record is None:
                status = MappingStatus.UNMAPPED.value
                entrez_id = None
                approved_symbol = gene.approved_symbol
                note = "Canonical gene ID not found in Ensembl 112 reference"
                unmapped_count += 1
            else:
                raw_entrez = ref_record.entrez_gene_id.strip() if ref_record.entrez_gene_id else ""
                approved_symbol = ref_record.gene_symbol or gene.approved_symbol

                if not raw_entrez:
                    status = MappingStatus.UNMAPPED.value
                    entrez_id = None
                    note = "No Entrez Gene ID associated in Ensembl 112 reference"
                    unmapped_count += 1
                elif "|" in raw_entrez:
                    # One-to-many mapping in reference
                    status = MappingStatus.ONE_TO_MANY.value
                    entrez_id = raw_entrez.split("|")[0].strip()
                    note = f"Multiple Entrez IDs found ({raw_entrez}); primary selected: {entrez_id}"
                    one_to_many_count += 1
                    mapped_count += 1
                    entrez_counts[entrez_id] = entrez_counts.get(entrez_id, 0) + 1
                else:
                    entrez_id = raw_entrez
                    mapped_count += 1
                    entrez_counts[entrez_id] = entrez_counts.get(entrez_id, 0) + 1
                    if not is_stat_valid:
                        status = MappingStatus.INVALID_STAT.value
                        note = f"Invalid Stage 7 statistical status: {gene.result_status} (stat={gene.wald_statistic})"
                        invalid_stat_count += 1
                    else:
                        status = MappingStatus.MAPPED.value
                        note = None

            intermediate.append((gene, entrez_id, status, note))

        # Second pass: flag duplicate Entrez targets
        duplicate_entrez_count = 0
        for gene, entrez_id, status, note in intermediate:
            final_status = status
            final_note = note
            if entrez_id and entrez_counts.get(entrez_id, 0) > 1:
                duplicate_entrez_count += 1
                if status == MappingStatus.MAPPED.value:
                    final_status = MappingStatus.DUPLICATE_ENTREZ.value
                    final_note = f"Multiple Ensembl genes map to Entrez {entrez_id}"

            records.append(
                PathwayGeneMappingRecord(
                    canonical_gene_id=gene.canonical_gene_id,
                    approved_symbol=gene.approved_symbol,
                    entrez_gene_id=entrez_id,
                    mapping_status=final_status,
                    source_reference_version=self.gene_reference.reference_version,
                    notes=final_note,
                )
            )

        return PathwayGeneMappingResult(
            mappings=tuple(records),
            total_input_genes=len(gene_results),
            mapped_count=mapped_count,
            unmapped_count=unmapped_count,
            duplicate_entrez_count=duplicate_entrez_count,
            one_to_many_count=one_to_many_count,
            invalid_stat_count=invalid_stat_count,
            reference_version=self.gene_reference.reference_version,
            reference_checksum=self.gene_reference.checksum,
        )
