"""Execution and path-translation utilities for WSL2 and Linux environments."""

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


def run_wsl(cmd: str) -> subprocess.CompletedProcess:
    """Run command in WSL dnt_bridge environment from either Windows or native Linux."""
    if sys.platform == "win32":
        full_cmd = f"source /home/chari/miniconda3/bin/activate dnt_bridge && {cmd}"
        return subprocess.run(
            ["wsl.exe", "-e", "bash", "-c", full_cmd],
            capture_output=True,
            text=True,
            check=False,
        )
    else:
        full_cmd = f"source /home/chari/miniconda3/bin/activate dnt_bridge 2>/dev/null || true; {cmd}"
        return subprocess.run(
            ["bash", "-c", full_cmd],
            capture_output=True,
            text=True,
            check=False,
        )
