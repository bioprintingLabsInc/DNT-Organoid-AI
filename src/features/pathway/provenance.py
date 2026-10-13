"""Provenance tracking and execution audit recording for Stage 8C."""

import subprocess
from typing import Any

from src.differential_expression.models import DifferentialExpressionContrastResult
from .models import PathwayAnalysisProvenance


def get_current_git_commit() -> str:
    """Return HEAD commit hash or fallback string."""
    try:
        res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if res.returncode == 0:
            return res.stdout.strip()
    except Exception:
        pass
    return "unknown"


def build_pathway_analysis_provenance(
    contrast: DifferentialExpressionContrastResult,
    ranking_metric: str = "wald_statistic",
    gene_reference_version: str = "Ensembl_112",
    pathway_databases: tuple[str, ...] = ("GO_BP", "KEGG"),
    software_packages: dict[str, str] | None = None,
    stage7_source_commit: str = "ca26d34e6c0263e743425cc7d6752d0e5773a9fe",
) -> PathwayAnalysisProvenance:
    """Construct full traceability provenance record for Stage 8C."""
    current_commit = get_current_git_commit()
    pkgs = software_packages or {
        "clusterProfiler": "4.14.6",
        "org.Hs.eg.db": "3.20.0",
        "fgsea": "1.32.4",
        "KEGGREST": "1.46.0",
        "R": "4.4.3",
        "Bioconductor": "3.20",
    }

    return PathwayAnalysisProvenance(
        study_id=contrast.study_id,
        experiment_id=contrast.experiment_id,
        contrast_id=contrast.contrast_id,
        agent_name=contrast.agent_name,
        agent_identifier=contrast.agent_identifier,
        concentration=contrast.concentration_or_dose,
        concentration_unit=contrast.concentration_or_dose_unit,
        exposure_duration=contrast.exposure_duration,
        duration_unit=contrast.exposure_duration_unit,
        developmental_stage=contrast.developmental_age_or_stage_at_exposure,
        treatment_sample_ids=contrast.treatment_sample_ids,
        control_sample_ids=contrast.control_sample_ids,
        gene_reference_version=gene_reference_version,
        ranking_metric=ranking_metric,
        pathway_databases=pathway_databases,
        pathway_database_version="Bioconductor 3.20 / org.Hs.eg.db 3.20.0 / KEGG.db",
        software_packages=pkgs,
        stage7_source_commit=stage7_source_commit,
        stage8c_source_commit=current_commit,
    )
