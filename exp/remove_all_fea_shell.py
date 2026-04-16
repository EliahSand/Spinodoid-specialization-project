#!/usr/bin/env python3
"""Remove all FEA_shell directories under the GNN raw samples tree.

Usage:
  python3 exp/remove_all_fea_shell.py
  python3 exp/remove_all_fea_shell.py --dry-run
  python3 exp/remove_all_fea_shell.py --root Matlab/GNN/data/raw/samples
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

DEFAULT_SAMPLES_ROOT = "Matlab/GNN/data/raw/samples"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Delete all FEA_shell folders under the GNN raw samples directory."
    )
    p.add_argument(
        "--root",
        default=DEFAULT_SAMPLES_ROOT,
        help="Root folder to scan recursively (default: %(default)s).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print folders that would be deleted without deleting them.",
    )
    return p.parse_args()


def resolve_root(root_arg: str) -> Path:
    root = Path(root_arg).expanduser()
    if root.is_absolute():
        return root
    repo_root = Path(__file__).resolve().parents[3]
    return (repo_root / root).resolve()


def find_fea_shell_dirs(scan_root: Path) -> list[Path]:
    if not scan_root.is_dir():
        return []
    return sorted(p for p in scan_root.rglob("FEA_shell") if p.is_dir())


def main() -> None:
    args = parse_args()
    scan_root = resolve_root(args.root)
    fea_dirs = find_fea_shell_dirs(scan_root)

    print("Scan root:", scan_root)
    print("Found %d FEA_shell director%s." % (len(fea_dirs), "y" if len(fea_dirs) == 1 else "ies"))

    if not fea_dirs:
        return

    removed = 0
    failed = 0
    for fea_dir in fea_dirs:
        if args.dry_run:
            print("[DRY ]", fea_dir)
            removed += 1
            continue

        try:
            shutil.rmtree(str(fea_dir))
            print("[DEL ]", fea_dir)
            removed += 1
        except Exception as exc:
            print("[FAIL] %s (%s)" % (fea_dir, exc))
            failed += 1

    mode = "dry-run" if args.dry_run else "delete"
    print("Done (%s): removed=%d, failed=%d." % (mode, removed, failed))


if __name__ == "__main__":
    main()
