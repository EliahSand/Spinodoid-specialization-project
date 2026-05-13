#!/usr/bin/env python
"""Batch run shell spinodal simulations in Abaqus (Python 2/3 compatible).

Usage:
    abaqus python exp/batch_run_spinodal_shell.py --roots Matlab/results [--overwrite] [--dry-run]
"""

from __future__ import print_function
import argparse
import codecs
import json
import os
import shutil
import subprocess
import time

HEAVY_EXTS = set([
    '.odb', '.cae', '.sim', '.prt', '.res', '.stt',
    '.msg', '.dat', '.com', '.log', '.lck', '.sta', '.jnl',
])


def parse_args():
    p = argparse.ArgumentParser(description='Batch run shell spinodal simulations.')
    p.add_argument('--roots', nargs='+', default=['Matlab/results'],
                   help='Root directories to search for mesh_manifest.json files.')
    p.add_argument('--abaqus-cmd', default='abaqus',
                   help='Abaqus command (default: abaqus)')
    p.add_argument('--script', default='exp/run_spinodal_shell_static.py',
                   help='Path to the shell runner script.')
    p.add_argument('--packaged-results-root', default=None,
                   help='Optional dataset/samples root. If a run already has '
                        'midpoint_results_shell.csv there, skip it unless --overwrite is used.')
    p.add_argument('--overwrite', action='store_true',
                   help='Re-run even if the target ODB already exists.')
    p.add_argument('--dry-run', action='store_true',
                   help='Print commands without executing.')
    p.add_argument('--max-parallel', type=int, default=1,
                   help='Number of Abaqus jobs to run concurrently (default: 1, i.e. serial).')
    p.add_argument('--poll-interval', type=float, default=2.0,
                   help='Seconds between polls of in-flight workers (default: 2.0).')
    return p.parse_args()


def find_manifests(roots):
    manifests = []
    for root in roots:
        root_path = os.path.expanduser(root)
        if not os.path.isdir(root_path):
            continue
        for dirpath, _dirnames, filenames in os.walk(root_path):
            if 'mesh_manifest.json' in filenames:
                manifests.append(os.path.join(dirpath, 'mesh_manifest.json'))
    return sorted(manifests)


def resolve_shell_inp_path(manifest_path, meta):
    mask_rel = meta.get('mask') or meta.get('mat')
    if not mask_rel:
        return None
    default_output = os.path.splitext(os.path.basename(mask_rel))[0] + '_shell.inp'
    manifest_output = meta.get('output_shell') or meta.get('output') or default_output
    return os.path.abspath(os.path.join(os.path.dirname(manifest_path), manifest_output))


def existing_odb_path(inp_path):
    job_name = os.path.splitext(os.path.basename(inp_path))[0] + '_shell_job'
    fea_dir = os.path.join(os.path.dirname(inp_path), 'FEA_shell')
    return os.path.join(fea_dir, job_name + '.odb')


def expected_midplane_csv_path(inp_path):
    fea_dir = os.path.join(os.path.dirname(inp_path), 'FEA_shell')
    return os.path.join(fea_dir, 'midplane_results_shell.csv')


def expected_packaged_midpoint_csv_path(packaged_results_root, manifest_path):
    if not packaged_results_root:
        return None
    run_name = os.path.basename(os.path.dirname(manifest_path))
    return os.path.join(os.path.abspath(packaged_results_root), run_name, 'midpoint_results_shell.csv')


def package_midplane_csv(packaged_results_root, manifest_path, mid_csv_path, overwrite=False, dry_run=False):
    if not packaged_results_root or not os.path.isfile(mid_csv_path):
        return False, None
    dst = expected_packaged_midpoint_csv_path(packaged_results_root, manifest_path)
    sample_dir = os.path.dirname(dst)
    if os.path.isfile(dst):
        if not overwrite:
            return False, dst
        if not dry_run:
            _safe_remove(dst)
    if dry_run:
        return True, dst
    if not os.path.isdir(sample_dir):
        os.makedirs(sample_dir)
    shutil.move(mid_csv_path, dst)
    return True, dst


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
            if ext == '.inp' and name.lower().endswith('_job.inp'):
                remove_file = True
            if remove_file:
                if not dry_run:
                    _safe_remove(path)
                removed += 1
        elif os.path.isdir(path) and name.lower().endswith('.simdir'):
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


def _post_process_finished_job(job, args):
    """Cleanup + package after a successful Abaqus run. Returns (cleaned, removed, cleanup_skip, packaged)."""
    mid_csv_path = job['mid_csv_path']
    fea_shell_dir = os.path.dirname(mid_csv_path)
    if not os.path.isfile(mid_csv_path):
        print('[WARN] %s: midplane_results_shell.csv not found; keeping Abaqus artifacts.' % job['inp_path'])
        return (0, 0, 1, 0)
    n = cleanup_fea_shell_dir(fea_shell_dir, dry_run=False)
    print('[CLEAN] %s: removed %d item(s).' % (fea_shell_dir, n))
    packaged = 0
    ok, dst = package_midplane_csv(
        args.packaged_results_root, job['manifest_path'], mid_csv_path,
        overwrite=args.overwrite, dry_run=False)
    if ok:
        packaged = 1
        print('[PACK] %s -> %s' % (mid_csv_path, dst))
    return (1, n, 0, packaged)


def _run_jobs_pool(jobs, args):
    """Run queued jobs with up to args.max_parallel concurrent Abaqus processes.

    Each job dict: {'manifest_path', 'inp_path', 'mid_csv_path'}.
    Returns aggregate (cleaned_runs, removed_items, cleanup_skips, packaged_runs).
    """
    cleaned_runs = 0
    removed_items = 0
    cleanup_skips = 0
    packaged_runs = 0

    in_flight = []  # list of dicts: {'proc', 'job', 't0'}
    next_idx = 0
    n_jobs = len(jobs)
    max_parallel = max(1, int(args.max_parallel))

    try:
        while next_idx < n_jobs or in_flight:
            # Launch new workers up to capacity.
            while len(in_flight) < max_parallel and next_idx < n_jobs:
                job = jobs[next_idx]
                cmd = '"%s" cae noGUI="%s" -- "%s"' % (args.abaqus_cmd, args.script, job['inp_path'])
                print('[RUN ] (%d/%d) %s' % (next_idx + 1, n_jobs, cmd))
                proc = subprocess.Popen(cmd, shell=True)
                in_flight.append({'proc': proc, 'job': job, 't0': time.time()})
                next_idx += 1

            # Poll for completion.
            time.sleep(args.poll_interval)
            still_running = []
            for entry in in_flight:
                ret = entry['proc'].poll()
                if ret is None:
                    still_running.append(entry)
                    continue
                dt = time.time() - entry['t0']
                job = entry['job']
                if ret != 0:
                    print('[FAIL] %s (exit %s, %.1fs)' % (job['inp_path'], ret, dt))
                else:
                    print('[DONE] %s (%.1fs)' % (job['inp_path'], dt))
                    c, r, s, p = _post_process_finished_job(job, args)
                    cleaned_runs += c
                    removed_items += r
                    cleanup_skips += s
                    packaged_runs += p
            in_flight = still_running
    except KeyboardInterrupt:
        print('\n[INTERRUPT] Terminating %d in-flight job(s)...' % len(in_flight))
        for entry in in_flight:
            try:
                entry['proc'].terminate()
            except Exception:
                pass
        # Best-effort wait so terminated children don't leave zombies.
        for entry in in_flight:
            try:
                entry['proc'].wait()
            except Exception:
                pass
        raise

    return cleaned_runs, removed_items, cleanup_skips, packaged_runs


def main():
    args = parse_args()
    manifests = find_manifests(args.roots)
    if not manifests:
        print('No mesh_manifest.json found under: %s' % ', '.join(args.roots))
        return

    batch_start = time.time()
    cleaned_runs = 0
    removed_items = 0
    cleanup_skips = 0
    packaged_runs = 0

        # --- Planning pass: decide skip / package-only / queue-for-run per manifest ---
    jobs_to_run = []

    for manifest_path in manifests:
        try:
            with codecs.open(manifest_path, 'r', encoding='utf-8-sig') as fh:
                meta = json.load(fh)
        except Exception as exc:
            print('[SKIP] %s: unable to read/parse JSON (%s).' % (manifest_path, exc))
            continue

        inp_path = resolve_shell_inp_path(manifest_path, meta)
        if inp_path is None:
            print('[SKIP] %s: missing mask/mat entry.' % manifest_path)
            continue
        if not os.path.isfile(inp_path):
            print('[SKIP] %s: shell INP not found.' % inp_path)
            continue

        packaged_csv_path = expected_packaged_midpoint_csv_path(args.packaged_results_root, manifest_path)
        if packaged_csv_path is not None and os.path.isfile(packaged_csv_path) and not args.overwrite:
            print('[SKIP] %s: packaged result exists (%s). Use --overwrite to rerun.' %
                  (inp_path, packaged_csv_path))
            continue

        odb_path = existing_odb_path(inp_path)
        mid_csv_path = expected_midplane_csv_path(inp_path)
        marker_path = None
        if os.path.isfile(mid_csv_path):
            marker_path = mid_csv_path
        elif os.path.isfile(odb_path):
            marker_path = odb_path
        if marker_path is not None and not args.overwrite:
            if marker_path == mid_csv_path and packaged_csv_path is not None:
                ok, dst = package_midplane_csv(
                    args.packaged_results_root, manifest_path, mid_csv_path,
                    overwrite=False, dry_run=args.dry_run)
                if ok:
                    packaged_runs += 1
                    print('[PACK] %s -> %s' % (mid_csv_path, dst))
                else:
                    print('[SKIP] %s: result exists (%s). Use --overwrite to rerun.' % (inp_path, marker_path))
                continue
            print('[SKIP] %s: result exists (%s). Use --overwrite to rerun.' % (inp_path, marker_path))
            continue
        
        jobs_to_run.append({
            'manifest_path': manifest_path,
            'inp_path': inp_path,
            'mid_csv_path': mid_csv_path,
        })

    print('Planned: %d job(s) to run, max_parallel=%d.' % (len(jobs_to_run), args.max_parallel))
    if args.dry_run:
        for j in jobs_to_run:
            print('[DRY ] would run: "%s" cae noGUI="%s" -- "%s"' %
                  (args.abaqus_cmd, args.script, j['inp_path']))
    elif jobs_to_run:
        c, r, s, p = _run_jobs_pool(jobs_to_run, args)
        cleaned_runs += c
        removed_items += r
        cleanup_skips += s
        packaged_runs += p

       

    batch_dt = time.time() - batch_start
    print('Batch complete in %.1f seconds (%.2f minutes)' % (batch_dt, batch_dt / 60.0))
    print('Immediate cleanup: cleaned %d run(s), removed %d item(s), skipped %d run(s), packaged %d run(s).' %
          (cleaned_runs, removed_items, cleanup_skips, packaged_runs))


if __name__ == '__main__':
    main()
