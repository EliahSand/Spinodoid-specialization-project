#!/usr/bin/env python3
"""Remove all solid .inp files under spesicifed directory

Usage:
  python3 cleaunUp/remove_all_sample_inp_solid.py
  python3 cleaunUp/remove_all_sample_inp_solid.py --dry-run
  python3 cleaunUp/remove_all_sample_inp_solid.py --root Matlab/defectPrediction/results
"""

from __future__ import annotations

import argparse
from pathlib import Path

DEFAULT_SAMPLES_ROOT = "Matlab/defectPrediction/results/"
TARGET_FILENAME = "sheet.inp"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Delete all solid inp files under the chosen directory."
    )
    parser.add_argument(
        "--root",
        default=DEFAULT_SAMPLES_ROOT,
        help="Root folder to scan recursively (default: %(default)s).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would be deleted without deleting them.",
    )
    return parser.parse_args()


def resolve_root(root_arg: str) -> Path:
    root = Path(root_arg).expanduser()
    if root.is_absolute():
        return root

    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return (parent / root).resolve()

    return (Path.cwd() / root).resolve()


def find_sample_mat_files(scan_root: Path) -> list[Path]:
    if not scan_root.is_dir():
        return []
    return sorted(path for path in scan_root.rglob(TARGET_FILENAME) if path.is_file())


def main() -> None:
    args = parse_args()
    scan_root = resolve_root(args.root)
    sample_mat_files = find_sample_mat_files(scan_root)

    print("Scan root:", scan_root)
    print(
        "Found %d %s file%s."
        % (
            len(sample_mat_files),
            TARGET_FILENAME,
            "" if len(sample_mat_files) == 1 else "s",
        )
    )

    if not sample_mat_files:
        return

    removed = 0
    failed = 0
    for sample_mat_file in sample_mat_files:
        if args.dry_run:
            print("[DRY ]", sample_mat_file)
            removed += 1
            continue

        try:
            sample_mat_file.unlink()
            print("[DEL ]", sample_mat_file)
            removed += 1
        except Exception as exc:
            print("[FAIL] %s (%s)" % (sample_mat_file, exc))
            failed += 1

    mode = "dry-run" if args.dry_run else "delete"
    print("Done (%s): removed=%d, failed=%d." % (mode, removed, failed))


if __name__ == "__main__":
    main()
