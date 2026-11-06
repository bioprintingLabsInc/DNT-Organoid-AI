"""Step 7B: Bulk Treatment-versus-Matched-Control Differential Expression v1."""

from .design import DesignResolutionResult, DesignResolver, matrix_rank
from .limma_runner import (
    execute_limma_trend_r,
    prepare_continuous_expression_matrix,
    run_limma_trend_contrast,
    run_limma_trend_pipeline,
)
from .models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    GeneDifferentialExpressionResult,
)
from .pipeline import run_contrast_differential_expression, run_differential_expression
from .r_bridge import (
    LOCKED_BIOC_VERSION,
    LOCKED_DESEQ2_VERSION,
    LOCKED_R_VERSION,
    execute_deseq2_contrast_r,
)

__all__ = [
    "DesignResolutionResult",
    "DesignResolver",
    "DifferentialExpressionContrastResult",
    "DifferentialExpressionDataset",
    "GeneDifferentialExpressionResult",
    "LOCKED_BIOC_VERSION",
    "LOCKED_DESEQ2_VERSION",
    "LOCKED_R_VERSION",
    "execute_deseq2_contrast_r",
    "execute_limma_trend_r",
    "matrix_rank",
    "prepare_continuous_expression_matrix",
    "run_contrast_differential_expression",
    "run_differential_expression",
    "run_limma_trend_contrast",
    "run_limma_trend_pipeline",
]
