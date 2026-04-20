#!/usr/bin/env python3
"""Run defect-only MAT -> INP conversion via the existing batch converter.

How to run from the repo root:
    python Matlab/defectPrediction/run_batch_defect_mat_to_inp.py

Run this after:
    run('Matlab/defectPrediction/run_batch_defect_generation.m')
"""

from pathlib import Path
import shlex
import subprocess
import sys


def main():
    repo_root = Path(__file__).resolve().parents[2]
    defect_root = repo_root / "Matlab" / "defectPrediction" / "results"
    manifests = sorted(defect_root.rglob("mesh_manifest.json"))
    if not manifests:
        raise SystemExit(
            "No generated defect manifests were found under Matlab/defectPrediction/results."
        )

    cmd = [
        sys.executable,
        "exp/batch_mat_to_inp.py",
        "--roots",
        "Matlab/defectPrediction/results",
    ]

    print(
        "Running solid MAT -> INP conversion for defect cases under "
        "Matlab/defectPrediction/results"
    )
    print("[CMD ] %s" % " ".join(shlex.quote(part) for part in cmd))

    result = subprocess.run(cmd, cwd=repo_root)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
