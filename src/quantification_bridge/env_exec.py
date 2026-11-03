"""Execution and path-translation utilities for WSL2 and Linux environments.

Enforces cross-laboratory portability:
1. No hardcoded user directories or machine-specific paths.
2. Dynamic conda / mamba environment discovery.
3. Native POSIX and WSL2 Windows path normalization.
4. Support for environment variable overrides (DNT_CONDA_ENV, DNT_CONDA_ACTIVATE).
"""

import os
import re
import subprocess
import sys
from pathlib import Path


def to_wsl_path(p: str | Path) -> str:
    """Convert any Windows or Linux Path object into a POSIX path suitable for WSL/Linux commands."""
    path_obj = Path(p).resolve()
    posix = path_obj.as_posix()
    if posix.startswith("/"):
        return posix
    match = re.match(r"^([a-zA-Z]):/(.*)", posix)
    if match:
        drive = match.group(1).lower()
        rest = match.group(2)
        return f"/mnt/{drive}/{rest}"
    return posix


def get_conda_activation_snippet(env_name: str | None = None) -> str:
    """Generate portable bash shell snippet to discover and activate the target conda environment."""
    target_env = env_name or os.environ.get("DNT_CONDA_ENV", "dnt_bridge")
    custom_activate = os.environ.get("DNT_CONDA_ACTIVATE")

    if custom_activate:
        return f"source {custom_activate} {target_env} 2>/dev/null || true; "

    # Portable multi-path discovery
    return (
        f'if [ -n "$CONDA_PREFIX" ] && [ -d "$CONDA_PREFIX/bin" ]; then '
        f'  PATH="$CONDA_PREFIX/bin:$PATH"; '
        f'elif command -v conda >/dev/null 2>&1; then '
        f'  eval "$(conda shell.bash hook 2>/dev/null)" && conda activate {target_env} 2>/dev/null || true; '
        f'else '
        f'  for c_base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" "$HOME/miniforge3" "/opt/conda" "/opt/miniconda3" "/usr/local/miniconda3"; do '
        f'    if [ -f "$c_base/bin/activate" ]; then '
        f'      source "$c_base/bin/activate" {target_env} 2>/dev/null && break || true; '
        f'    fi; '
        f'  done; '
        f'fi; '
    )


def run_wsl(cmd: str, env_name: str | None = None) -> subprocess.CompletedProcess:
    """Run command in the target bioinformatics environment via WSL2 (on Windows) or bash (on Linux)."""
    activate_prefix = get_conda_activation_snippet(env_name)
    full_cmd = f"{activate_prefix}{cmd}"

    if sys.platform == "win32":
        return subprocess.run(
            ["wsl.exe", "-e", "bash", "-c", full_cmd],
            capture_output=True,
            text=True,
            check=False,
        )
    else:
        return subprocess.run(
            ["bash", "-c", full_cmd],
            capture_output=True,
            text=True,
            check=False,
        )
