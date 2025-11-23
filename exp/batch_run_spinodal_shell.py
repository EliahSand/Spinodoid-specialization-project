#!/usr/bin/env python3
"""Batch run shell spinodal simulations using run_spinodal_shell_static.py."""

import argparse
import json
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional


def parse_args():
    p = argparse.ArgumentParser(description='Batch run shell spinodal simulations.')
    p.add_argument('--roots', nargs='+', default=['Matlab/results'],
                   help='Root directories to search for mesh_manifest.json files.')
    p.add_argument('--abaqus-cmd', default='abaqus',
                   help='Abaqus command (default: abaqus)')
    p.add_argument('--script', default='exp/run_spinodal_shell_static.py',
                   help='Path to the shell runner script.')
    p.add_argument('--overwrite', action='store_true',
                   help='Re-run even if the target ODB already exists.')
    p.add_argument('--dry-run', action='store_true',
                   help='Print commands without executing.')
    return p.parse_args()


def find_manifests(roots: Iterable[str]) -> List[Path]:
    manifests: List[Path] = []
    for root in roots:
        root_path = Path(root).expanduser()
        if not root_path.exists():
            continue
        manifests.extend(root_path.rglob('mesh_manifest.json'))
    return sorted(manifests)


def resolve_shell_inp_path(manifest_path: Path, meta: dict) -> Optional[Path]:
    mask_rel = meta.get('mask') or meta.get('mat')
    if not mask_rel:
        return None
    default_output = f"{Path(mask_rel).stem}_shell.inp"
    manifest_output = meta.get('output_shell') or meta.get('output') or default_output
    return (manifest_path.parent / manifest_output).resolve()


def existing_odb_path(inp_path: Path) -> Path:
    job_name = inp_path.stem + '_shell_job'
    fea_dir = inp_path.parent / 'FEA_shell'
    return fea_dir / f'{job_name}.odb'


def main():
    args = parse_args()
    manifests = find_manifests(args.roots)
    if not manifests:
        print('No mesh_manifest.json found under:', ', '.join(args.roots))
        return

    for manifest_path in manifests:
        try:
            meta = json.loads(manifest_path.read_text())
        except Exception:
            print(f'[SKIP] {manifest_path}: unable to read/parse JSON.')
            continue
        inp_path = resolve_shell_inp_path(manifest_path, meta)
        if inp_path is None:
            print(f'[SKIP] {manifest_path}: missing mask/mat entry.')
            continue
        if not inp_path.exists():
            print(f'[SKIP] {inp_path}: shell INP not found.')
            continue
        odb_path = existing_odb_path(inp_path)
        if odb_path.exists() and not args.overwrite:
            print(f'[SKIP] {inp_path}: ODB exists ({odb_path}). Use --overwrite to rerun.')
            continue

        cmd = [
            args.abaqus_cmd,
            'cae',
            f'noGUI={args.script}',
            '--',
            str(inp_path),
        ]
        print('[RUN ]', ' '.join(cmd))
        if args.dry_run:
            continue
        ret = subprocess.call(cmd)
        if ret != 0:
            print(f'[FAIL] {inp_path} (exit {ret})')
        else:
            print(f'[DONE] {inp_path}')


if __name__ == '__main__':
    main()
