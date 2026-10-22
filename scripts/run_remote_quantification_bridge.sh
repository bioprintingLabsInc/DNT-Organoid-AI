#!/usr/bin/env bash
# ==============================================================================
# Governed Remote Execution Wrapper for PRJNA802072 Quantification Bridge
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "=== PRJNA802072 Quantification Bridge Remote Runner ==="
echo "Repository root: ${REPO_ROOT}"

# Activate conda environment
ENV_NAME="dnt_bridge"
if [[ -f "${HOME}/miniconda3/bin/activate" ]]; then
    source "${HOME}/miniconda3/bin/activate" "${ENV_NAME}"
elif [[ -f "${HOME}/anaconda3/bin/activate" ]]; then
    source "${HOME}/anaconda3/bin/activate" "${ENV_NAME}"
elif command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}"
else
    echo "WARNING: Conda activation script not found in standard paths. Relying on active PATH."
fi

# Verify core dependencies
echo "Verifying bioinformatics toolchain on node..."
which fastp || { echo "ERROR: fastp missing from PATH"; exit 1; }
which STAR || { echo "ERROR: STAR missing from PATH"; exit 1; }
which featureCounts || { echo "ERROR: featureCounts missing from PATH"; exit 1; }
which samtools || { echo "ERROR: samtools missing from PATH"; exit 1; }
which infer_experiment.py || { echo "ERROR: infer_experiment.py missing from PATH"; exit 1; }

# Parse arguments
THREADS="${1:-16}"

echo "Starting pipeline with ${THREADS} threads..."
cd "${REPO_ROOT}"
python3 "${REPO_ROOT}/scripts/run_remote_quantification_bridge.py" --threads "${THREADS}" "$@"

echo "=== Remote execution step complete ==="
