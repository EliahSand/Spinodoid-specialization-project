#!/usr/bin/env python3
"""Quarantine duplicate raw GNN samples by (thickness, angle, rng_seed).

Default mode is a dry-run. Use --apply to move duplicate folders under
Matlab/GNN/data/raw/quarantine_duplicate_samples while preserving the same
relative trXX/angYYY/run-folder layout.
"""

from __future__ import annotations

import argparse
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_SAMPLES_ROOT = "Matlab/GNN/data/raw/samples"
DEFAULT_DATASET_ROOT = "Matlab/GNN/data/dataset/samples"
DEFAULT_HYBRID_ROOT = "Matlab/GNN/data/dataset_hybrid/samples"
DEFAULT_QUARANTINE_ROOT = "Matlab/GNN/data/raw/quarantine_duplicate_samples"

SEED_RE = re.compile(r"rng_seed:\s*(\d+)")


@dataclass(frozen=True)
class RunInfo:
    key: Tuple[str, str, int]
    path: Path
    rel: Path
    run_name: str
    in_dataset: bool
    in_hybrid: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Keep one raw run per (trXX, angYYY, rng_seed) and quarantine duplicates."
    )
    parser.add_argument(
        "--samples-root",
        default=DEFAULT_SAMPLES_ROOT,
        help="Raw samples root to curate.",
    )
    parser.add_argument(
        "--dataset-root",
        default=DEFAULT_DATASET_ROOT,
        help="Packaged Abaqus dataset samples root used to prefer canonical folders.",
    )
    parser.add_argument(
        "--hybrid-root",
        default=DEFAULT_HYBRID_ROOT,
        help="Hybrid dataset samples root used as a secondary canonical preference.",
    )
    parser.add_argument(
        "--quarantine-root",
        default=DEFAULT_QUARANTINE_ROOT,
        help="Destination root for duplicate folders when --apply is used.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Move duplicate folders. Omit for dry-run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Explicit dry-run flag for readability. This is the default unless --apply is used.",
    )
    parser.add_argument(
        "--examples",
        type=int,
        default=12,
        help="Number of duplicate actions to print.",
    )
    return parser.parse_args()


def read_seed(log_path: Path) -> Optional[int]:
    try:
        text = log_path.read_text(errors="replace")
    except OSError:
        return None
    match = SEED_RE.search(text)
    if not match:
        return None
    return int(match.group(1))


def sample_names(root: Path) -> set[str]:
    if not root.is_dir():
        return set()
    return {path.name for path in root.iterdir() if path.is_dir()}


def iter_run_dirs(samples_root: Path) -> Iterable[Path]:
    for tr_dir in sorted(samples_root.glob("tr*")):
        if not tr_dir.is_dir():
            continue
        for ang_dir in sorted(tr_dir.glob("ang*")):
            if not ang_dir.is_dir():
                continue
            for run_dir in sorted(ang_dir.iterdir()):
                if run_dir.is_dir():
                    yield run_dir


def discover_runs(
    samples_root: Path,
    dataset_root: Path,
    hybrid_root: Path,
) -> tuple[list[RunInfo], list[tuple[Path, str]]]:
    packaged = sample_names(dataset_root)
    hybrid = sample_names(hybrid_root)

    runs: list[RunInfo] = []
    invalid: list[tuple[Path, str]] = []
    for run_dir in iter_run_dirs(samples_root):
        rel = run_dir.relative_to(samples_root)
        parts = rel.parts
        if len(parts) != 3:
            invalid.append((rel, "unexpected_depth"))
            continue
        tr_label, ang_label, run_name = parts
        log_path = run_dir / "run_log.txt"
        if not log_path.is_file():
            invalid.append((rel, "missing_run_log"))
            continue
        seed = read_seed(log_path)
        if seed is None:
            invalid.append((rel, "missing_rng_seed"))
            continue
        runs.append(
            RunInfo(
                key=(tr_label, ang_label, seed),
                path=run_dir,
                rel=rel,
                run_name=run_name,
                in_dataset=run_name in packaged,
                in_hybrid=run_name in hybrid,
            )
        )
    return runs, invalid


def canonical_sort_key(run: RunInfo) -> tuple[int, int, str]:
    return (
        0 if run.in_dataset else 1,
        0 if run.in_hybrid else 1,
        str(run.rel),
    )


def choose_canonical(runs: Sequence[RunInfo]) -> tuple[RunInfo, list[RunInfo]]:
    ordered = sorted(runs, key=canonical_sort_key)
    return ordered[0], ordered[1:]


def build_actions(runs: Sequence[RunInfo]) -> tuple[list[RunInfo], list[RunInfo], Dict[Tuple[str, str, int], list[RunInfo]]]:
    grouped: Dict[Tuple[str, str, int], list[RunInfo]] = defaultdict(list)
    for run in runs:
        grouped[run.key].append(run)

    keep: list[RunInfo] = []
    quarantine: list[RunInfo] = []
    duplicate_groups: Dict[Tuple[str, str, int], list[RunInfo]] = {}
    for key, group in grouped.items():
        if len(group) == 1:
            keep.append(group[0])
            continue
        canonical, duplicates = choose_canonical(group)
        keep.append(canonical)
        quarantine.extend(duplicates)
        duplicate_groups[key] = group
    return keep, quarantine, duplicate_groups


def print_summary(
    runs: Sequence[RunInfo],
    invalid: Sequence[tuple[Path, str]],
    keep: Sequence[RunInfo],
    quarantine: Sequence[RunInfo],
    duplicate_groups: Dict[Tuple[str, str, int], list[RunInfo]],
    examples: int,
) -> None:
    folders_per_combo = Counter((run.key[0], run.key[1]) for run in runs)
    seeds_per_combo: dict[tuple[str, str], set[int]] = defaultdict(set)
    for run in runs:
        seeds_per_combo[(run.key[0], run.key[1])].add(run.key[2])

    duplicate_hist = Counter(len(group) for group in duplicate_groups.values())
    unique_seed_hist = Counter(len(seeds) for seeds in seeds_per_combo.values())
    folder_hist = Counter(folders_per_combo.values())

    print("Raw duplicate curation summary")
    print("  parsed_runs: %d" % len(runs))
    print("  invalid_runs: %d" % len(invalid))
    print("  unique_tr_angle_seed: %d" % len(keep))
    print("  duplicate_groups: %d" % len(duplicate_groups))
    print("  duplicate_extra_folders: %d" % len(quarantine))
    print("  keep_already_packaged: %d" % sum(run.in_dataset for run in keep))
    print("  keep_already_hybrid: %d" % sum(run.in_hybrid for run in keep))
    print("  quarantine_already_packaged: %d" % sum(run.in_dataset for run in quarantine))
    print("  quarantine_already_hybrid: %d" % sum(run.in_hybrid for run in quarantine))
    print("  duplicate_multiplicity_hist: %s" % format_counter(duplicate_hist))
    print("  unique_seeds_per_combo_hist: %s" % format_counter(unique_seed_hist))
    print("  folders_per_combo_hist: %s" % format_counter(folder_hist))

    if invalid:
        print("Invalid run examples:")
        for rel, reason in invalid[:examples]:
            print("  %s: %s" % (rel, reason))

    if quarantine:
        print("Quarantine examples:")
        for run in sorted(quarantine, key=lambda item: str(item.rel))[:examples]:
            print(
                "  %s seed=%d packaged=%s hybrid=%s"
                % (run.rel, run.key[2], run.in_dataset, run.in_hybrid)
            )


def format_counter(counter: Counter) -> str:
    if not counter:
        return "(none)"
    return ", ".join("%s:%s" % (key, counter[key]) for key in sorted(counter))


def apply_quarantine(samples_root: Path, quarantine_root: Path, duplicates: Sequence[RunInfo]) -> None:
    for run in sorted(duplicates, key=lambda item: str(item.rel)):
        destination = quarantine_root / run.rel
        if destination.exists():
            raise FileExistsError(
                "Quarantine destination already exists: %s. Move it aside before rerunning --apply."
                % destination
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(samples_root / run.rel), str(destination))


def main() -> None:
    args = parse_args()
    if args.apply and args.dry_run:
        raise SystemExit("Choose only one of --apply or --dry-run.")

    samples_root = Path(args.samples_root).expanduser()
    dataset_root = Path(args.dataset_root).expanduser()
    hybrid_root = Path(args.hybrid_root).expanduser()
    quarantine_root = Path(args.quarantine_root).expanduser()

    if not samples_root.is_dir():
        raise SystemExit("Samples root not found: %s" % samples_root)

    runs, invalid = discover_runs(samples_root, dataset_root, hybrid_root)
    keep, quarantine, duplicate_groups = build_actions(runs)
    print_summary(runs, invalid, keep, quarantine, duplicate_groups, args.examples)

    if not args.apply:
        print("Dry-run only. Re-run with --apply to move duplicate folders.")
        return

    apply_quarantine(samples_root, quarantine_root, quarantine)
    print("Moved %d duplicate folder(s) to %s" % (len(quarantine), quarantine_root))


if __name__ == "__main__":
    main()
