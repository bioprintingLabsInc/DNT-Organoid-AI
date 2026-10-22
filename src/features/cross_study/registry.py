"""Registry constructor and query services for Stage 8A.1 Condition-Level Gene Response.

Governing Invariants:
1. One record per gene × contrast (no aggregation, no collapsing).
2. Exact numerical preservation of Stage 7 statistics.
3. Multiple doses/timepoints from the same study do NOT increase independent study count.
4. Single-study evidence assigned CrossStudyStatus.SINGLE_STUDY_ONLY.
5. Zero meta-analysis or cross-study pooling executed in Stage 8A.1.
6. Zero DNT outcome labels consumed.
"""

from collections import defaultdict
from collections.abc import Sequence
from datetime import datetime, timezone
from typing import Any

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
)
from src.features.cross_study.builder import (
    build_records_for_contrast,
    normalize_independent_study_id,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
)


def build_condition_gene_response_registry(
    contrasts: Sequence[DifferentialExpressionContrastResult] | DifferentialExpressionDataset,
    registry_id: str = "reg_cross_study_gene_response_v1",
    registry_version: str = "1.0.0",
    study_family_map: dict[str, str] | None = None,
    model_type_map: dict[str, str] | None = None,
) -> ConditionGeneResponseRegistry:
    """Build a complete, deterministic ConditionGeneResponseRegistry from Stage 7 contrast results.

    Parameters
    ----------
    contrasts : Sequence[DifferentialExpressionContrastResult] | DifferentialExpressionDataset
        One or more Stage 7 differential expression contrast results.
    registry_id : str
        Identifier for the registry instance.
    registry_version : str
        Governed version of the registry format.
    study_family_map : dict[str, str] | None
        Optional mapping from study_id to study_family_id for multi-experiment series.
    model_type_map : dict[str, str] | None
        Optional mapping from contrast_id to model/organoid system name.
    """
    contrast_list: list[DifferentialExpressionContrastResult] = []
    if isinstance(contrasts, DifferentialExpressionDataset):
        contrast_list = list(contrasts.contrast_results)
    else:
        contrast_list = list(contrasts)

    if not contrast_list:
        return ConditionGeneResponseRegistry(
            records=(),
            registry_id=registry_id,
            registry_version=registry_version,
            created_at=datetime.now(timezone.utc).isoformat(),
            provenance={"total_contrasts": 0, "total_records": 0},
        )

    family_map = study_family_map or {}
    model_map = model_type_map or {}

    # 1. First pass: Determine independent study count per compound
    compound_to_studies: dict[str, set[str]] = defaultdict(set)
    for c in contrast_list:
        comp_id = c.agent_identifier or c.agent_name or "UNKNOWN_COMPOUND"
        fam_id = family_map.get(c.study_id)
        ind_study_id = normalize_independent_study_id(c.study_id, fam_id)
        compound_to_studies[comp_id].add(ind_study_id)

    # 2. Second pass: Construct ConditionGeneResponseRecords with accurate study status
    all_records: list[ConditionGeneResponseRecord] = []
    for c in contrast_list:
        comp_id = c.agent_identifier or c.agent_name or "UNKNOWN_COMPOUND"
        study_count = len(compound_to_studies[comp_id])

        if study_count <= 1:
            cross_status = CrossStudyStatus.SINGLE_STUDY_ONLY
        else:
            # Future Stage 8A.2 will evaluate comparability; never auto-promote
            cross_status = CrossStudyStatus.MULTI_STUDY_CONTEXT_UNDER_REVIEW

        fam_id = family_map.get(c.study_id)
        mod_type = model_map.get(c.contrast_id)

        contrast_records = build_records_for_contrast(
            contrast=c,
            study_family_id=fam_id,
            model_type=mod_type,
            cross_study_status=cross_status,
            registry_version=registry_version,
        )
        all_records.extend(contrast_records)

    # 3. Deterministic sorting: by study_id, contrast_id, canonical_ensembl_gene_id
    all_records.sort(
        key=lambda r: (
            r.study_id,
            r.contrast_id,
            r.canonical_ensembl_gene_id,
        )
    )

    provenance_info: dict[str, Any] = {
        "registry_id": registry_id,
        "registry_version": registry_version,
        "total_contrasts": len(contrast_list),
        "total_records": len(all_records),
        "unique_canonical_genes_count": len({r.canonical_ensembl_gene_id for r in all_records}),
        "unique_contrasts": [c.contrast_id for c in contrast_list],
        "unique_compounds": list(compound_to_studies.keys()),
        "independent_study_counts_by_compound": {
            comp: len(studies) for comp, studies in compound_to_studies.items()
        },
        "meta_analysis_executed": False,
        "dnt_labels_used": False,
    }

    return ConditionGeneResponseRegistry(
        records=tuple(all_records),
        registry_id=registry_id,
        registry_version=registry_version,
        created_at=datetime.now(timezone.utc).isoformat(),
        provenance=provenance_info,
    )
