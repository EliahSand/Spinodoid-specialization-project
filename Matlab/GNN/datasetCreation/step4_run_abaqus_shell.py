#!/usr/bin/env python
"""Step 4 wrapper: run Abaqus shell batch and package midpoint CSVs.

How to run:
  abaqus cae noGUI=Matlab/GNN/datasetCreation/step4_run_abaqus_shell.py --
  abaqus cae noGUI=Matlab/GNN/datasetCreation/step4_run_abaqus_shell.py -- --dry-run
 

Writes:
  Matlab/GNN/data/dataset/samples/<run_name>/midpoint_results_shell.csv
"""

from __future__ import print_function

import argparse
import os
from shutil import move
import sys
import time


PY2 = sys.version_info[0] == 2
if PY2:
    text_type = unicode  # noqa: F821 (only defined in Python 2)
else:
    text_type = str

DEFAULT_SAMPLES_ROOT = "Matlab/GNN/data/raw/samples"
DEFAULT_DATASET_ROOT = "Matlab/GNN/data/dataset/samples"
DEFAULT_ABAQUS_CMD = "abaqus"

HEAVY_EXTS = set([
    ".odb", ".cae", ".sim", ".prt", ".res", ".stt",
    ".msg", ".dat", ".com", ".log", ".lck", ".sta", ".jnl",
])


def _extract_user_cli_args(argv):
    """Allow launch via both python and `abaqus cae noGUI=... -- ...`."""
    if "--" in argv:
        idx = argv.index("--")
        return argv[idx + 1:]
    return argv


def parse_args():
    parser = argparse.ArgumentParser(description="Run Abaqus shell batch for GNN dataset.")
    parser.add_argument("--dry-run", action="store_true", help="Print command only.")
    parser.add_argument(
        "--package-only",
        action="store_true",
        help="Recovery mode: skip Abaqus and package existing midplane_results_shell.csv files.",
    )
    argv = _extract_user_cli_args(sys.argv[1:])
    args, _unknown = parser.parse_known_args(argv)
    return args


def _safe_remove(path):
    try:
        os.remove(path)
    except OSError:
        pass


def cleanup_fea_shell_dir(fea_shell_dir, dry_run):
    if not os.path.isdir(fea_shell_dir):
        return 0
    removed = 0
    for name in os.listdir(fea_shell_dir):
        path = os.path.join(fea_shell_dir, name)
        if os.path.isfile(path):
            ext = os.path.splitext(name)[1].lower()
            remove_file = ext in HEAVY_EXTS
            if ext == ".inp" and name.lower().endswith("_job.inp"):
                remove_file = True
            if remove_file:
                if not dry_run:
                    _safe_remove(path)
                removed += 1
        elif os.path.isdir(path) and name.lower().endswith(".simdir"):
            if dry_run:
                removed += 1
            else:
                for root, dirs, files in os.walk(path, topdown=False):
                    for fn in files:
                        _safe_remove(os.path.join(root, fn))
                    for dn in dirs:
                        dpath = os.path.join(root, dn)
                        try:
                            os.rmdir(dpath)
                        except OSError:
                            pass
                try:
                    os.rmdir(path)
                except OSError:
                    pass
                removed += 1
    return removed


def _to_text(value):
    if value is None:
        return ""
    if isinstance(value, text_type):
        return value
    if PY2:
        try:
            return value.decode("utf-8")
        except Exception:
            return text_type(value)
    try:
        if isinstance(value, bytes):
            return value.decode("utf-8")
    except Exception:
        pass
    return str(value)


def package_processed_sample(run_name, csv_path, dataset_root, dry_run):
    sample_dir = os.path.join(dataset_root, run_name)
    midpoint_dst = os.path.join(sample_dir, "midpoint_results_shell.csv")
    src_csv = csv_path
    if not os.path.isfile(src_csv):
        if os.path.isfile(midpoint_dst):
            src_csv = midpoint_dst
        else:
            return False, "missing_midpoint_csv"

    if dry_run:
        return True, "dry_run"

    if not os.path.isdir(sample_dir):
        os.makedirs(sample_dir)
    if os.path.abspath(src_csv) != os.path.abspath(midpoint_dst):
        if os.path.isfile(midpoint_dst):
            _safe_remove(midpoint_dst)
        move(src_csv, midpoint_dst)
    return True, "ok"


def iter_midplane_csvs(samples_root):
    for dirpath, _dirnames, filenames in os.walk(samples_root):
        if "midplane_results_shell.csv" in filenames:
            yield os.path.join(dirpath, "midplane_results_shell.csv")


def infer_run_name_from_midpoint_csv(csv_path):
    parent = os.path.dirname(csv_path)
    if os.path.basename(parent).lower() == "fea_shell":
        return os.path.basename(os.path.dirname(parent))
    return os.path.basename(parent)


def resolve_script_path():
    # Normal Python execution.
    if "__file__" in globals():
        return os.path.abspath(__file__)

    # Abaqus/CAE noGUI often executes via execfile without __file__.
    for arg in sys.argv:
        txt = _to_text(arg).strip().strip('"').strip("'")
        low = txt.lower()
        if "nogui=" in low:
            rhs = txt.split("=", 1)[1].strip().strip('"').strip("'")
            if os.path.isfile(rhs):
                return os.path.abspath(rhs)

    # Fallbacks.
    if sys.argv and os.path.isfile(sys.argv[0]):
        return os.path.abspath(sys.argv[0])

    guess = os.path.abspath(os.path.join(os.getcwd(), "Matlab", "GNN", "pipeline", "step4_run_abaqus_shell.py"))
    if os.path.isfile(guess):
        return guess

    return os.path.abspath(os.getcwd())


def resolve_repo_root(script_path):
    if os.path.isfile(script_path):
        candidate = os.path.abspath(os.path.join(os.path.dirname(script_path), "..", "..", ".."))
        batch_script = os.path.join(candidate, "exp", "batch_run_spinodal_shell.py")
        if os.path.isfile(batch_script):
            return candidate

    cwd_candidate = os.path.abspath(os.getcwd())
    batch_script = os.path.join(cwd_candidate, "exp", "batch_run_spinodal_shell.py")
    if os.path.isfile(batch_script):
        return cwd_candidate

    raise RuntimeError(
        "Could not resolve project root containing exp/batch_run_spinodal_shell.py. "
        "Run from repo root or keep default project layout."
    )


def run_batch_shell(batch_script_path, roots, abaqus_cmd, packaged_results_root, dry_run):
    argv_backup = list(sys.argv)
    try:
        sys.argv = [
            batch_script_path,
            "--roots",
            roots,
            "--abaqus-cmd",
            abaqus_cmd,
            "--packaged-results-root",
            packaged_results_root,
        ]
        if dry_run:
            sys.argv.append("--dry-run")
        namespace = {
            "__name__": "__main__",
            "__file__": batch_script_path,
        }
        if PY2:
            execfile(batch_script_path, namespace)
        else:
            with open(batch_script_path, "rb") as fh:
                code = compile(fh.read(), batch_script_path, "exec")
            exec(code, namespace)
    finally:
        sys.argv = argv_backup


def main():
    args = parse_args()
    script_path = resolve_script_path()
    repo_root = resolve_repo_root(script_path)
    script = os.path.join(repo_root, "exp", "batch_run_spinodal_shell.py")
    dataset_root = os.path.abspath(os.path.join(repo_root, DEFAULT_DATASET_ROOT))
    samples_root = os.path.abspath(os.path.join(repo_root, DEFAULT_SAMPLES_ROOT))

    batch_time = time.time()

    print(
        "Running: %s python %s --roots %s --abaqus-cmd %s --packaged-results-root %s%s"
        % (
            DEFAULT_ABAQUS_CMD,
            script,
            samples_root,
            DEFAULT_ABAQUS_CMD,
            dataset_root,
            " --dry-run" if args.dry_run else "",
        )
    )
    if not args.package_only:
        run_batch_shell(
            batch_script_path=script,
            roots=samples_root,
            abaqus_cmd=DEFAULT_ABAQUS_CMD,
            packaged_results_root=dataset_root,
            dry_run=args.dry_run,
        )

        cleaned_runs = 0
        removed_items = 0
        for csv_path in iter_midplane_csvs(samples_root):
            fea_shell_dir = os.path.dirname(csv_path)
            n = cleanup_fea_shell_dir(fea_shell_dir, dry_run=args.dry_run)
            if n > 0:
                cleaned_runs += 1
                removed_items += n
        print("Post-run cleanup: cleaned %d run(s), removed %d item(s)." % (cleaned_runs, removed_items))
    else:
        print("Package-only mode: skipping Abaqus run and cleanup.")

    package_csv_paths = list(iter_midplane_csvs(samples_root))
    if args.package_only and not package_csv_paths:
        package_csv_paths = list(iter_midplane_csvs(dataset_root))

    packaged = 0
    skipped = 0
    for csv_path in package_csv_paths:
        run_name = infer_run_name_from_midpoint_csv(csv_path)
        ok, reason = package_processed_sample(
            run_name=run_name,
            csv_path=csv_path,
            dataset_root=dataset_root,
            dry_run=args.dry_run,
        )
        if ok:
            packaged += 1
        else:
            skipped += 1
            print("[SKIP] %s: %s" % (run_name, reason))
    print(
        "Dataset packaging: packaged=%d, skipped=%d, output_root=%s"
        % (packaged, skipped, dataset_root)
    )

    final_time = time.time() - batch_time
    print('Batch complete in %.1f seconds and %.2f minutes' % (final_time, final_time / 60.0))


# This wrapper is only used as an entry script. Abaqus/CAE noGUI executes it
# via execfile(..., __main__.__dict__), so rely on a one-shot sentinel rather
# than fragile __name__ / argv detection.
if not globals().get("_STEP4_RUN_ABAQUS_SHELL_ALREADY_EXECUTED", False):
    globals()["_STEP4_RUN_ABAQUS_SHELL_ALREADY_EXECUTED"] = True
    main()
