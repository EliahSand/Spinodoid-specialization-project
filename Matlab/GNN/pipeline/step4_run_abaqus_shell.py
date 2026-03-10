#!/usr/bin/env python
"""Step 4 wrapper: run Abaqus shell batch and package midpoint targets.

How to run:
  abaqus cae noGUI=Matlab/GNN/pipeline/step4_run_abaqus_shell.py --
  abaqus cae noGUI=Matlab/GNN/pipeline/step4_run_abaqus_shell.py -- --dry-run
 

Writes:
  Matlab/GNN/data/dataset/samples/<run_name>/midpoint_results_shell.csv
  Matlab/GNN/data/dataset/samples/<run_name>/target_midpoint_curvature.json
  Matlab/GNN/data/dataset/samples/<run_name>/target_midpoint_curvature.csv
  Matlab/GNN/data/dataset/samples/<run_name>/meta.json
"""

from __future__ import print_function

import argparse
import csv
import json
import math
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
        help="Skip Abaqus run and package targets from existing midplane_results_shell.csv files.",
    )
    argv = _extract_user_cli_args(sys.argv[1:])
    args, _unknown = parser.parse_known_args(argv)
    return args


def _is_finite(v):
    try:
        return math.isfinite(v)
    except AttributeError:
        return not (math.isnan(v) or math.isinf(v))


def _finite_or_none(v):
    if v is None:
        return None
    try:
        vf = float(v)
    except Exception:
        return None
    if _is_finite(vf):
        return vf
    return None


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


def _normalize_header(h):
    txt = _to_text(h)
    return txt.lstrip(u"\ufeff").strip()


def _open_csv_read(path):
    if PY2:
        return open(path, "rb")
    return open(path, "r", newline="", encoding="utf-8-sig")


def _read_midplane_rows(csv_path):
    rows_out = []
    with _open_csv_read(csv_path) as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return rows_out
        keys = [_normalize_header(h) for h in header]
        for row in reader:
            if not row:
                continue
            if len(row) < len(keys):
                row = list(row) + [""] * (len(keys) - len(row))
            out = {}
            for i, key in enumerate(keys):
                out[key] = _to_text(row[i]) if i < len(row) else ""
            rows_out.append(out)
    return rows_out


def _load_midplane_deformed_profile(csv_path):
    y_def = []
    z_def = []
    for row in _read_midplane_rows(csv_path):
        try:
            y = float(row["Y"])
            z = float(row["Z"])
            u2 = float(row["U2"])
            u3 = float(row["U3"])
        except Exception:
            continue
        yd = y + u2
        zd = z + u3
        if _is_finite(yd) and _is_finite(zd):
            y_def.append(yd)
            z_def.append(zd)
    return y_def, z_def


def _kasa_circle_fit(y_vals, z_vals):
    """Algebraic circle fit (Kasa), same model as MATLAB check_curvature_single_run."""
    try:
        import numpy as np
    except Exception:
        return None

    y = np.asarray(y_vals, dtype=float)
    z = np.asarray(z_vals, dtype=float)
    if y.size < 3:
        return None

    a = np.column_stack((2.0 * y, 2.0 * z, np.ones_like(y)))
    b = y ** 2 + z ** 2
    try:
        try:
            res = np.linalg.lstsq(a, b, rcond=None)
        except TypeError:
            res = np.linalg.lstsq(a, b)
        x = res[0]
    except Exception:
        return None

    yc = float(x[0])
    zc = float(x[1])
    c0 = float(x[2])
    r2 = c0 + yc ** 2 + zc ** 2
    if (not _is_finite(r2)) or r2 <= 0:
        return None

    radius = math.sqrt(r2)
    if radius <= 0:
        return None
    rr = np.hypot(y - yc, z - zc) - radius
    rmse = float(math.sqrt(float(np.mean(rr ** 2)))) if rr.size else None
    return {
        "kappa_circle_1_per_m": 1.0 / radius,
        "radius_circle_m": radius,
        "circle_center_y_m": yc,
        "circle_center_z_m": zc,
        "circle_fit_rmse_m": rmse,
    }


def compute_midpoint_curvature(csv_path):
    y, z = _load_midplane_deformed_profile(csv_path)
    n = len(y)
    out = {
        "n_points": n,
        "kappa_circle_1_per_m": None,
        "radius_circle_m": None,
        "circle_center_y_m": None,
        "circle_center_z_m": None,
        "circle_fit_rmse_m": None,
        "kappa_sagitta_1_per_m": None,
        "radius_sagitta_m": None,
        "delta_sagitta_m": None,
        "chord_length_m": None,
    }
    if n < 3:
        return out

    pts = sorted(zip(y, z), key=lambda t: t[0])
    y = [p[0] for p in pts]
    z = [p[1] for p in pts]
    y0, z0 = y[0], z[0]
    y1, z1 = y[-1], z[-1]
    dy = y1 - y0
    dz = z1 - z0
    chord = max(y) - min(y)
    out["chord_length_m"] = _finite_or_none(chord)
    vnorm = math.hypot(dy, dz)
    if chord > 0 and vnorm > 0:
        signed = [((yy - y0) * dz - (zz - z0) * dy) / vnorm for yy, zz in zip(y, z)]
        delta = max(abs(v) for v in signed) if signed else None
        if delta is not None and delta > 0:
            radius_s = (chord ** 2) / (8.0 * delta) + delta / 2.0
            out["delta_sagitta_m"] = _finite_or_none(delta)
            out["radius_sagitta_m"] = _finite_or_none(radius_s)
            out["kappa_sagitta_1_per_m"] = _finite_or_none(1.0 / radius_s if radius_s > 0 else None)

    circle = _kasa_circle_fit(y, z)
    if circle:
        out.update(circle)
    return out


def write_curvature_csv(curvature_dict, out_path):
    if PY2:
        fh = open(out_path, "wb")
    else:
        fh = open(out_path, "w", newline="", encoding="utf-8")
    try:
        writer = csv.writer(fh)
        headers = [
            "n_points",
            "kappa_circle_1_per_m",
            "radius_circle_m",
            "circle_center_y_m",
            "circle_center_z_m",
            "circle_fit_rmse_m",
            "kappa_sagitta_1_per_m",
            "radius_sagitta_m",
            "delta_sagitta_m",
            "chord_length_m",
        ]
        writer.writerow(headers)
        writer.writerow([curvature_dict.get(k) for k in headers])
    finally:
        fh.close()


def package_processed_sample(run_name, csv_path, dataset_root, dry_run):
    sample_dir = os.path.join(dataset_root, run_name)
    midpoint_dst = os.path.join(sample_dir, "midpoint_results_shell.csv")
    src_csv = csv_path
    if not os.path.isfile(src_csv):
        if os.path.isfile(midpoint_dst):
            src_csv = midpoint_dst
        else:
            return False, "missing_midpoint_csv"

    curvature = compute_midpoint_curvature(src_csv)
    target_dst = os.path.join(sample_dir, "target_midpoint_curvature.json")
    target_csv_dst = os.path.join(sample_dir, "target_midpoint_curvature.csv")
    meta_dst = os.path.join(sample_dir, "meta.json")
    sample_mat_path = os.path.join(sample_dir, "sample.mat")
    sample_mat_exists = os.path.isfile(sample_mat_path)

    if dry_run:
        return True, "dry_run"

    if not os.path.isdir(sample_dir):
        os.makedirs(sample_dir)
    if os.path.abspath(src_csv) != os.path.abspath(midpoint_dst):
        if os.path.isfile(midpoint_dst):
            _safe_remove(midpoint_dst)
        move(src_csv, midpoint_dst)
    src_csv = midpoint_dst

    with open(target_dst, "w") as fh:
        json.dump(curvature, fh, indent=2)
    write_curvature_csv(curvature, target_csv_dst)
    with open(meta_dst, "w") as fh:
        json.dump({
            "run_name": run_name,
            "source_midplane_csv": src_csv,
            "midpoint_results_shell_csv": midpoint_dst,
            "sample_mat_path": sample_mat_path,
            "sample_mat_exists": sample_mat_exists,
            "target_name": "kappa_circle_1_per_m",
            "target_curvature_json": target_dst,
            "target_curvature_csv": target_csv_dst,
        }, fh, indent=2)
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


def _load_batch_module(batch_script_path):
    if PY2:
        import imp
        return imp.load_source("batch_run_spinodal_shell", batch_script_path)

    import importlib.machinery
    loader = importlib.machinery.SourceFileLoader("batch_run_spinodal_shell", batch_script_path)
    return loader.load_module()


def run_batch_shell(batch_script_path, roots, abaqus_cmd, dry_run):
    mod = _load_batch_module(batch_script_path)
    if not hasattr(mod, "main"):
        raise RuntimeError("batch_run_spinodal_shell.py has no main()")

    argv_backup = list(sys.argv)
    try:
        sys.argv = [batch_script_path, "--roots", roots, "--abaqus-cmd", abaqus_cmd]
        if dry_run:
            sys.argv.append("--dry-run")
        mod.main()
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
        "Running: %s python %s --roots %s --abaqus-cmd %s%s"
        % (
            DEFAULT_ABAQUS_CMD,
            script,
            DEFAULT_SAMPLES_ROOT,
            DEFAULT_ABAQUS_CMD,
            " --dry-run" if args.dry_run else "",
        )
    )
    if not args.package_only:
        run_batch_shell(
            batch_script_path=script,
            roots=DEFAULT_SAMPLES_ROOT,
            abaqus_cmd=DEFAULT_ABAQUS_CMD,
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
    print('Batch complete in %.1f seconds and %.2f minutes' % (batch_time, batch_time/60))


if __name__ == "__main__":
    main()
