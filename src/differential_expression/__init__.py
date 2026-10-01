"""Step 7B: Bulk Treatment-versus-Matched-Control Differential Expression v1."""

from .design import DesignResolutionResult, DesignResolver, matrix_rank
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
    "matrix_rank",
    "run_contrast_differential_expression",
    "run_differential_expression",
]
