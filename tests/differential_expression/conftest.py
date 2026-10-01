"""Test fixtures and synthetic data generators for Differential Expression v1."""

from typing import Any
import random
import pytest

from src.gene_harmonization.models import (
    CanonicalCollision,
    GeneIdentifierType,
    GeneMappingStatus,
    HarmonizationSummary,
    HarmonizedDataset,
    HarmonizedGene,
    Status as HarmonizationStatus,
)
from src.normalization.models import (
    NormalizationCohort,
    NormalizedCohortResult,
    NormalizedCountMatrix,
    NormalizedDataset,
    ResponseEligibility,
    SampleSizeFactor,
    Status as NormalizationStatus,
)
from src.qc.matrix import CountMatrix
from src.response_builder.models import (
    ContrastExposure,
    ContrastStatus,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Status as ContrastDatasetStatus,
)


def create_synthetic_de_fixtures(
    n_treatment: int = 3,
    n_control: int = 3,
    n_genes: int = 25,
    custom_gene_counts: dict[str, tuple[list[int], list[int]]] | None = None,
    size_factor_val: float = 1.0,
    treatment_bio_reps: list[str] | None = None,
    control_bio_reps: list[str] | None = None,
    biological_sources: list[str] | None = None,
    treatment_source: str = "cell_line_A",
    control_source: str = "cell_line_A",
    sample_sources: dict[str, str] | None = None,
    contrast_status: ContrastStatus = ContrastStatus.ELIGIBLE,
    developmental_age: str = "Day 30",
    control_developmental_age: str = "Day 30",
    exposures: list[ContrastExposure] | None = None,
    include_unmapped_gene: bool = False,
    include_colliding_gene: bool = False,
    use_diagnostic_as_raw: bool = False,
    invalid_size_factor: float | None = None,
) -> tuple[MolecularResponseContrast, NormalizedDataset, HarmonizedDataset, dict[str, Any]]:
    """Generate complete synthetic fixtures for testing Step 7B differential expression."""
    trt_sids = [f"trt_s{i+1}" for i in range(n_treatment)]
    ctrl_sids = [f"ctrl_s{i+1}" for i in range(n_control)]
    all_sids = ctrl_sids + trt_sids

    trt_reps = treatment_bio_reps or [f"rep_trt_{i+1}" for i in range(n_treatment)]
    ctrl_reps = control_bio_reps or [f"rep_ctrl_{i+1}" for i in range(n_control)]

    # Gene definitions
    gene_ids: list[str] = []
    harmonized_genes: list[HarmonizedGene] = []
    collisions: list[CanonicalCollision] = []

    # Standard background genes
    for i in range(1, n_genes + 1):
        can_id = f"ENSG{i:011d}"
        orig_id = f"GENE_{i}"
        sym = f"SYM{i}"
        gene_ids.append(orig_id)
        harmonized_genes.append(
            HarmonizedGene(
                source_index=i - 1,
                original_gene_id=orig_id,
                identifier_type=GeneIdentifierType.GENE_SYMBOL,
                canonical_gene_id=can_id,
                approved_symbol=sym,
                mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
                mapping_reason=None,
                candidate_canonical_ids=(can_id,),
                reference_id="Ensembl",
                reference_version="112",
                source_asset_id="asset_01",
            )
        )

    if include_unmapped_gene:
        idx = len(gene_ids)
        gene_ids.append("UNMAPPED_GENE")
        harmonized_genes.append(
            HarmonizedGene(
                source_index=idx,
                original_gene_id="UNMAPPED_GENE",
                identifier_type=GeneIdentifierType.CUSTOM,
                canonical_gene_id=None,
                approved_symbol=None,
                mapping_status=GeneMappingStatus.UNMAPPED,
                mapping_reason="Not found in reference",
                candidate_canonical_ids=(),
                reference_id="Ensembl",
                reference_version="112",
                source_asset_id="asset_01",
            )
        )

    if include_colliding_gene:
        idx = len(gene_ids)
        gene_ids.append("COLLIDING_GENE_1")
        gene_ids.append("COLLIDING_GENE_2")
        can_coll = "ENSG_COLLISION_01"
        harmonized_genes.append(
            HarmonizedGene(
                source_index=idx,
                original_gene_id="COLLIDING_GENE_1",
                identifier_type=GeneIdentifierType.CUSTOM,
                canonical_gene_id=can_coll,
                approved_symbol="COLL1",
                mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
                mapping_reason=None,
                candidate_canonical_ids=(can_coll,),
                reference_id="Ensembl",
                reference_version="112",
                source_asset_id="asset_01",
            )
        )
        harmonized_genes.append(
            HarmonizedGene(
                source_index=idx + 1,
                original_gene_id="COLLIDING_GENE_2",
                identifier_type=GeneIdentifierType.CUSTOM,
                canonical_gene_id=can_coll,
                approved_symbol="COLL1",
                mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
                mapping_reason=None,
                candidate_canonical_ids=(can_coll,),
                reference_id="Ensembl",
                reference_version="112",
                source_asset_id="asset_01",
            )
        )
        collisions.append(
            CanonicalCollision(
                canonical_gene_id=can_coll,
                original_gene_ids=("COLLIDING_GENE_1", "COLLIDING_GENE_2"),
                source_indices=(idx, idx + 1),
                approved_symbol="COLL1",
            )
        )

    # Build count columns
    # Background counts: ~ 100 + i * 5
    counts_by_sample: dict[str, list[int]] = {s: [] for s in all_sids}
    for i, gid in enumerate(gene_ids):
        if custom_gene_counts and gid in custom_gene_counts:
            ctrl_vals, trt_vals = custom_gene_counts[gid]
            for s_idx, s in enumerate(ctrl_sids):
                val = ctrl_vals[s_idx % len(ctrl_vals)]
                counts_by_sample[s].append(val)
            for s_idx, s in enumerate(trt_sids):
                val = trt_vals[s_idx % len(trt_vals)]
                counts_by_sample[s].append(val)
        else:
            mean_levels = [50, 80, 120, 200, 350, 500, 800, 1200]
            mu = mean_levels[i % len(mean_levels)]
            gene_rng = random.Random(1000 + i)
            for s in all_sids:
                lam = gene_rng.gammavariate(10, mu / 10)
                val = max(1, int(round(gene_rng.gauss(lam, (mu / 10) ** 0.5))))
                counts_by_sample[s].append(val)

    raw_columns = tuple(tuple(counts_by_sample[s]) for s in all_sids)
    raw_matrix = CountMatrix(gene_ids=tuple(gene_ids), sample_ids=tuple(all_sids), columns=raw_columns)

    if use_diagnostic_as_raw:
        float_cols = tuple(tuple(float(v) for v in col) for col in raw_columns)
        raw_matrix = NormalizedCountMatrix(gene_ids=tuple(gene_ids), sample_ids=tuple(all_sids), columns=float_cols)

    # Size factors
    size_factors = []
    for sid in all_sids:
        sf = size_factor_val
        if invalid_size_factor is not None and sid == trt_sids[0]:
            sf = invalid_size_factor
        size_factors.append(
            SampleSizeFactor(
                sample_id=sid,
                size_factor=sf,
                raw_library_size=sum(counts_by_sample[sid]),
                normalized_library_size=float(sum(counts_by_sample[sid])),
                cohort_id="cohort_01",
                response_eligibility=ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL if sid in trt_sids else ResponseEligibility.CONTROL_BASELINE,
                matched_control_sample_ids=tuple(ctrl_sids) if sid in trt_sids else (),
                treatment_control_status="treatment" if sid in trt_sids else "control",
            )
        )

    cohort = NormalizationCohort(
        cohort_id="cohort_01",
        study_id="study_01",
        experiment_id="exp_01",
        sample_ids=tuple(all_sids),
        assay_ids=("assay_01",),
        rna_seq_modality="bulk_rna_seq",
        library_strategy="polyA",
        sequencing_platform="Illumina",
        strandedness="reverse",
        read_layout="paired_end",
        treatment_sample_ids=tuple(trt_sids),
        control_sample_ids=tuple(ctrl_sids),
        treatment_control_relationship_ids=("rel_01",),
    )

    cohort_result = NormalizedCohortResult(
        status=NormalizationStatus.PASS,
        cohort=cohort,
        findings=(),
        size_factors=tuple(size_factors),
        metrics=None,
        diagnostic_normalized_matrix=None,
        raw_count_matrix=raw_matrix,
        method="ratio",
        r_environment_info={"r_version": "4.4.3", "bioc_version": "3.20", "deseq2_version": "1.46.0"},
        execution_provenance={},
    )

    normalized_dataset = NormalizedDataset(
        status=NormalizationStatus.PASS,
        findings=(),
        cohort_results=(cohort_result,),
        source_asset_id="asset_01",
        reference_identity="Ensembl_112",
        config_version="1.0.0",
        normalization_version="1.0.0",
    )

    harmonization_summary = HarmonizationSummary(
        total_genes=len(gene_ids),
        uniquely_mapped_genes=sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED),
        unmapped_genes=sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.UNMAPPED),
        ambiguous_genes=sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.AMBIGUOUS),
        invalid_genes=0,
        canonical_collision_genes=len(collisions),
        colliding_source_gene_count=sum(len(c.original_gene_ids) for c in collisions),
        percentage_uniquely_mapped=100.0,
        percentage_unresolved=0.0,
    )

    harmonized_dataset = HarmonizedDataset(
        status=HarmonizationStatus.PASS,
        findings=(),
        summary=harmonization_summary,
        genes=tuple(harmonized_genes),
        collisions=tuple(collisions),
        matrix=raw_matrix if isinstance(raw_matrix, CountMatrix) else None,
        source_asset_id="asset_01",
        assay_id="assay_01",
        sample_ids=tuple(all_sids),
        reference_identity="Ensembl_112",
        reference_version="112",
        reference_checksum="dummy_checksum",
        config_version="1.0.0",
        config_checksum="dummy_config_checksum",
        harmonization_version="1.0.0",
    )

    # Sample metadata
    sample_meta: dict[str, Any] = {}
    for sid in all_sids:
        is_trt = sid in trt_sids
        if sample_sources and sid in sample_sources:
            src = sample_sources[sid]
        else:
            src = treatment_source if is_trt else control_source
        b_rep = trt_reps[trt_sids.index(sid)] if is_trt else ctrl_reps[ctrl_sids.index(sid)]
        age = developmental_age if is_trt else control_developmental_age
        sample_meta[sid] = {
            "sample_id": sid,
            "condition_id": "treat_01" if is_trt else "ctrl_01",
            "biological_replicate_id": b_rep,
            "biological_source_id": src,
            "developmental_age": age,
        }

    contrast_exposures = tuple(exposures or [
        ContrastExposure(
            exposure_id="exp_01",
            agent_name="Valproic acid",
            agent_identifier="CHEBI:30953",
            concentration_or_dose=200.0,
            concentration_or_dose_unit="uM",
            exposure_duration=24.0,
            exposure_duration_unit="hours",
            developmental_age_or_stage_at_exposure=developmental_age,
        )
    ])

    contrast = MolecularResponseContrast(
        contrast_id="contrast_001",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="treat_01",
        matched_control_condition_ids=("ctrl_01",),
        treatment_sample_ids=tuple(trt_sids),
        control_sample_ids=tuple(ctrl_sids),
        treatment_biological_replicate_ids=tuple(trt_reps),
        control_biological_replicate_ids=tuple(ctrl_reps),
        treatment_exposures=contrast_exposures,
        agent_name=contrast_exposures[0].agent_name if len(contrast_exposures) == 1 else None,
        agent_identifier=contrast_exposures[0].agent_identifier if len(contrast_exposures) == 1 else None,
        concentration_or_dose=contrast_exposures[0].concentration_or_dose if len(contrast_exposures) == 1 else None,
        concentration_or_dose_unit=contrast_exposures[0].concentration_or_dose_unit if len(contrast_exposures) == 1 else None,
        exposure_duration=contrast_exposures[0].exposure_duration if len(contrast_exposures) == 1 else None,
        exposure_duration_unit=contrast_exposures[0].exposure_duration_unit if len(contrast_exposures) == 1 else None,
        developmental_age_or_stage_at_exposure=developmental_age,
        treatment_collection_age_or_stage=developmental_age,
        control_collection_age_or_stage=control_developmental_age,
        treatment_collection_age_normalized=developmental_age,
        control_collection_age_normalized=control_developmental_age,
        biological_source_ids=tuple(biological_sources or [treatment_source, control_source]),
        status=contrast_status,
    )

    return contrast, normalized_dataset, harmonized_dataset, sample_meta
