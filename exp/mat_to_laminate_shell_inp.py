#!/usr/bin/env python3
"""Convert laminate MAT mask to a shell Abaqus INP.

Usage:
    python exp/mat_to_laminate_shell_inp.py --mat Matlab/results/.../sheet.mat
"""

import argparse

from mat_to_shell_inp import convert_mat_to_shell


def parse_args():
    p = argparse.ArgumentParser(description='Laminate MAT -> shell INP converter')
    p.add_argument('--mat', required=True, help='Path to .mat file containing sheetMask')
    p.add_argument('--var', default='sheetMask', help='Mask variable name (default: sheetMask)')
    p.add_argument('--spacing', type=float, default=None,
                   help='Optional voxel size override (meters). Usually omitted.')
    p.add_argument('--output', default=None, help='Output .inp path (default: mat stem + _shell.inp)')
    p.add_argument('--material', default='SPINODAL', help='Material name in INP')
    p.add_argument('--density', type=float, default=None, help='Optional density')
    p.add_argument('--elastic', type=float, nargs=2, metavar=('E', 'NU'), default=None,
                   help='Optional elastic pair (E, nu)')
    p.add_argument('--origin', type=float, nargs=3, default=None,
                   help='Optional XYZ origin override')
    p.add_argument('--peel', type=int, default=0, help='Optional peel voxels from all faces')
    p.add_argument('--base-thickness', type=float, default=None, help='Optional base shell thickness override')
    p.add_argument('--spin-thickness', type=float, default=None, help='Optional spin shell thickness override')
    return p.parse_args()


def main():
    args = parse_args()
    out_path, base_count, spin_count = convert_mat_to_shell(
        mat_path=args.mat,
        varname=args.var,
        spacing=args.spacing,
        origin=tuple(args.origin) if args.origin is not None else None,
        material=args.material,
        density=args.density,
        elastic=tuple(args.elastic) if args.elastic is not None else None,
        output_path=args.output,
        peel=args.peel,
        base_thickness=args.base_thickness,
        spin_thickness=args.spin_thickness,
    )
    print('Wrote laminate shell mesh to {} (BASE elements: {}, SPIN elements: {})'.format(
        out_path, base_count, spin_count))


if __name__ == '__main__':
    main()
