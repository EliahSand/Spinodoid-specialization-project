#!/usr/bin/env python3
"""Step 2 wrapper: reuse existing batch MAT->shell INP converter.

How to run:
  python3 Matlab/GNN/pipeline/step2_mat_to_shell_inp.py
  python3 Matlab/GNN/pipeline/step2_mat_to_shell_inp.py --dry-run
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

DEFAULT_SAMPLES_ROOT = "Matlab/GNN/data/raw/samples"
DEFAULT_PEEL = 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert MAT to shell INP for GNN dataset.")
    p.add_argument("--dry-run", action="store_true", help="Print command only.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[3]
    script = repo_root / "exp" / "batch_mat_to_shell_inp.py"

    # Reuse the same interpreter used to launch this wrapper (works on Windows without `python` alias).
    cmd = [
        sys.executable,
        str(script),
        "--roots",
        DEFAULT_SAMPLES_ROOT,
        "--peel",
        str(DEFAULT_PEEL),
    ]
    if args.dry_run:
        cmd.append("--dry-run")

    print("Running:", " ".join(cmd))
    ret = subprocess.call(cmd, cwd=str(repo_root))
    if ret != 0:
        raise SystemExit(ret)


if __name__ == "__main__":
    main()
