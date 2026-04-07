#!/usr/bin/env python
"""Batch run solid spinodal simulations in Abaqus (Python 2/3 compatible).

Usage:
    abaqus python exp/batch_run_spinodal_solid.py --roots Matlab/results [--overwrite] [--dry-run]
"""


from __future__ import print_function
import argparse
import codecs
import json
import os
import subprocess
import time


def parse_args():
    p = argparse.ArgumentParser(description='Batch run solid spinodal simulations.')
    p.add_argument('--roots', nargs='+', default=['Matlab/results'],
                   help='Root directories to search for mesh_manifest.json files.')
    p.add_argument('--abaqus-cmd', default='abaqus',
                   help='Abaqus command (default: abaqus)')
    p.add_argument('--script', default='exp/run_spinodal_static.py',
                   help='Path to the solid runner script.')
    p.add_argument('--overwrite', action='store_true',
                   help='Re-run even if the target ODB already exists.')
    p.add_argument('--dry-run', action='store_true',
                   help='Print commands without executing.')
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


def resolve_solid_inp_path(manifest_path, meta):
    mask_rel = meta.get('mask') or meta.get('mat')
    if not mask_rel:
        return None
    default_output = os.path.splitext(os.path.basename(mask_rel))[0] + '.inp'
    manifest_output = meta.get('output') or default_output
    return os.path.abspath(os.path.join(os.path.dirname(manifest_path), manifest_output))


def existing_odb_path(inp_path):
    job_name = os.path.splitext(os.path.basename(inp_path))[0] + '_job'
    fea_dir = os.path.join(os.path.dirname(inp_path), 'FEA')
    return os.path.join(fea_dir, job_name + '.odb')


def main():
    args = parse_args()
    manifests = find_manifests(args.roots)
    if not manifests:
        print('No mesh_manifest.json found under: %s' % ', '.join(args.roots))
        return

    batch_start = time.time()

    for manifest_path in manifests:
        try:
            with codecs.open(manifest_path, 'r', encoding='utf-8-sig') as fh:
                meta = json.load(fh)
        except Exception as exc:
            print('[SKIP] %s: unable to read/parse JSON (%s).' % (manifest_path, exc))
            continue

        inp_path = resolve_solid_inp_path(manifest_path, meta)
        if inp_path is None:
            print('[SKIP] %s: missing mask/mat entry.' % manifest_path)
            continue
        if not os.path.isfile(inp_path):
            print('[SKIP] %s: solid INP not found.' % inp_path)
            continue

        odb_path = existing_odb_path(inp_path)
        if os.path.isfile(odb_path) and not args.overwrite:
            print('[SKIP] %s: ODB exists (%s). Use --overwrite to rerun.' % (inp_path, odb_path))
            continue

        cmd = '"%s" cae noGUI="%s" -- "%s"' % (args.abaqus_cmd, args.script, inp_path)
        print('[RUN ] %s' % cmd)
        if args.dry_run:
            continue
        t0 = time.time()
        ret = subprocess.call(cmd, shell=True)
        dt = time.time() - t0
        if ret != 0:
            print('[FAIL] %s (exit %s, %.1fs)' % (inp_path, ret, dt))
        else:
            print('[DONE] %s (%.1fs)' % (inp_path, dt))

    batch_dt = time.time() - batch_start
    print('Batch complete in %.1f seconds (%.2f minutes)' % (batch_dt, batch_dt / 60.0))


if __name__ == '__main__':
    main()
