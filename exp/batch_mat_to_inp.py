#!/usr/bin/env python3
"""Batch converter that walks run folders and invokes mat_to_inp for each manifest."""

"""     how to run:
        python exp/batch_mat_to_inp.py --roots Matlab/results   
"""


import argparse
import json
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

from mat_to_inp import convert_mat_to_inp


def parse_args():
    parser = argparse.ArgumentParser(
        description='Scan manifest files and run mat_to_inp on each entry.')
    parser.add_argument('--roots', nargs='+', default=['Matlab/results'],
                        help='Root directories to search for mesh_manifest.json files.')
    parser.add_argument('--output-root', default=None,
                        help='Optional folder where generated INP files are mirrored. '
                             'Defaults to each run folder.')
    parser.add_argument('--material', default=None,
                        help='Override material name for all runs.')
    parser.add_argument('--density', type=float, default=None,
                        help='Override density value for all runs.')
    parser.add_argument('--elastic', type=float, nargs=2, metavar=('E', 'NU'), default=None,
                        help='Override elastic values (E, nu) for all runs.')
    parser.add_argument('--peel', type=int, default=0,
                        help='Peel value forwarded to mat_to_inp (default 0).')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing .inp files.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print the actions without writing files.')
    return parser.parse_args()


def find_manifests(roots: Iterable[str]) -> List[Path]:
    manifests: List[Path] = []
    for root in roots:
        root_path = Path(root).expanduser()
        if not root_path.exists():
            continue
        manifests.extend(root_path.rglob('mesh_manifest.json'))
    return sorted(manifests)


def resolve_relpath(manifest: Path, roots: Iterable[str]) -> Optional[Tuple[Path, Path]]:
    for root in roots:
        root_path = Path(root).expanduser().resolve()
        try:
            rel = manifest.parent.resolve().relative_to(root_path)
            return root_path, rel
        except ValueError:
            continue
    return None


def load_manifest(path: Path) -> dict:
    with path.open('r') as fh:
        return json.load(fh)


def main():
    args = parse_args()
    manifests = find_manifests(args.roots)
    if not manifests:
        print('No mesh_manifest.json files found under:', ', '.join(args.roots))
        return

    for manifest_path in manifests:
        manifest_dir = manifest_path.parent
        meta = load_manifest(manifest_path)
        mask_rel = meta.get('mask') or meta.get('mat')
        if not mask_rel:
            print(f'[SKIP] {manifest_path}: manifest missing "mask" entry.')
            continue
        mat_path = (manifest_dir / mask_rel).resolve()
        if not mat_path.exists():
            print(f'[SKIP] {manifest_path}: mask file {mat_path} missing.')
            continue
        var_name = meta.get('var')
        spacing = meta.get('spacing')
        origin = meta.get('origin')
        material = args.material or meta.get('material') or 'SPINODAL'
        density = args.density if args.density is not None else meta.get('density')
        elastic = args.elastic if args.elastic is not None else meta.get('elastic')

        default_output = Path(mask_rel).with_suffix('.inp').name
        manifest_output = meta.get('output')
        out_name = manifest_output or default_output

        if args.output_root:
            resolved = resolve_relpath(manifest_path, args.roots)
            if resolved:
                _, rel = resolved
                out_path = Path(args.output_root).expanduser().resolve() / rel / out_name
            else:
                out_path = Path(args.output_root).expanduser().resolve() / out_name
        else:
            out_path = manifest_dir / out_name

        if out_path.exists() and not args.overwrite:
            print(f'[SKIP] {out_path} already exists (use --overwrite to regenerate).')
            continue

        print(f'[RUN ] {mat_path} -> {out_path}')
        if args.dry_run:
            continue
        out_path.parent.mkdir(parents=True, exist_ok=True)
        convert_mat_to_inp(
            mat_path=str(mat_path),
            varname=var_name,
            spacing=spacing,
            origin=origin,
            material=material,
            density=density,
            elastic=tuple(elastic) if elastic is not None else None,
            output_path=str(out_path),
            peel=args.peel,
        )


if __name__ == '__main__':
    main()
