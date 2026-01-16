#!/usr/bin/env python3
"""Convert a MATLAB spinodal sheet mask into a shell element INP.

usage:

python exp/mat_to_shell_inp.py --mat Matlab/results/.../sheet.mat --var sheetMask --spacing <voxel_size_m>
"""

import argparse
import os
from typing import Optional, Tuple

import numpy as np

from mat_to_inp import (
    _detect_sheet_layers,
    _normalize_spacing,
    _resolve_origin,
    _resolve_spacing,
    load_mask_from_mat,
)


def parse_args():
    p = argparse.ArgumentParser(description='MATLAB mask -> shell INP converter')
    p.add_argument('--mat', required=True, help='Path to .mat file containing the mask')
    p.add_argument('--var', default=None,
                   help='Variable name inside the .mat file (default: auto-detect first 3D logical array)')
    p.add_argument('--spacing', type=float, default=None,
                   help='Voxel edge length (meters). If omitted, attempts to read "voxelSpacing" from the .mat')
    p.add_argument('--output', default=None, help='Output .inp path (default: mat stem + _shell.inp)')
    p.add_argument('--material', default='SPINODAL', help='Material name for *Material card')
    p.add_argument('--density', type=float, default=None, help='Optional density value')
    p.add_argument('--elastic', type=float, nargs=2, metavar=('E', 'NU'), default=None,
                   help='Optional elastic pair (E, nu)')
    p.add_argument('--base-thickness', type=float, default=None,
                   help='Override base shell thickness (meters)')
    p.add_argument('--spin-thickness', type=float, default=None,
                   help='Override spinodal shell thickness (meters, not added to base)')
    p.add_argument('--origin', type=float, nargs=3, default=None,
                   help='XYZ offset applied to node coordinates (default: from .mat or [0,0,0])')
    p.add_argument('--peel', type=int, default=0,
                   help='Number of voxels to strip from each face before meshing (default 0)')
    return p.parse_args()


def _write_id_lines(fh, ids, header, per_line=16):
    fh.write(header + '\n')
    line = []
    count = 0
    for item in ids:
        line.append(str(item))
        count += 1
        if count >= per_line:
            fh.write(', '.join(line) + '\n')
            line = []
            count = 0
    if line:
        fh.write(', '.join(line) + '\n')


def _infer_boundary_from_thickness(mat_data) -> Optional[int]:
    if 'zVoxelThickness' not in mat_data:
        return None
    zth = np.array(mat_data['zVoxelThickness']).astype(float).ravel()
    if zth.size < 2:
        return None
    first = zth[0]
    tol = max(1e-12, 1e-6 * max(1.0, abs(first)))
    diff_idx = np.where(np.abs(zth - first) > tol)[0]
    if diff_idx.size == 0:
        return None
    boundary = int(diff_idx[0])
    if boundary <= 0 or boundary >= zth.size:
        return None
    return boundary


def compute_thicknesses(mat_data, mask, spacing, layer_split) -> Tuple[float, float, int]:
    sx, sy, sz = _normalize_spacing(spacing)
    boundary = None
    if layer_split and layer_split.get('axis', None) == 0:
        boundary = int(layer_split.get('boundary', 0))
    if boundary is None:
        boundary = _infer_boundary_from_thickness(mat_data)

    zth = None
    if 'zVoxelThickness' in mat_data:
        zth = np.array(mat_data['zVoxelThickness']).astype(float).ravel()

    if boundary is None:
        # Fallback: first slice whose occupancy drops indicates spinodal start
        zdim = mask.shape[0]
        flat = mask.reshape(zdim, -1)
        full = flat.shape[1]
        threshold = 0.99 * full
        for idx, val in enumerate(flat.sum(axis=1)):
            if val < threshold:
                boundary = idx
                break
    if boundary is None:
        boundary = mask.shape[0] // 2
    boundary = max(0, min(int(boundary), mask.shape[0]))

    base_th = 0.0
    spin_th = 0.0
    if zth is not None and boundary <= zth.size:
        base_th = float(np.sum(zth[:boundary])) if boundary > 0 else 0.0
        spin_th = float(np.sum(zth[boundary:])) if boundary < zth.size else 0.0
    else:
        base_th = boundary * sz
        spin_th = max(mask.shape[0] - boundary, 0) * sz
    return base_th, spin_th, boundary


def write_shell_inp(path: str,
                    mask: np.ndarray,
                    spacing,
                    origin: Tuple[float, float, float],
                    material: str,
                    density: Optional[float] = None,
                    elastic: Optional[Tuple[float, float]] = None,
                    mat_data: Optional[dict] = None,
                    layer_split: Optional[dict] = None,
                    base_thickness: Optional[float] = None,
                    spin_thickness: Optional[float] = None) -> Tuple[int, int]:
    zdim, ydim, xdim = mask.shape
    ox, oy, oz = origin
    sx, sy, sz = _normalize_spacing(spacing)

    mat_for_thickness = mat_data if mat_data is not None else {}
    auto_base, auto_spin, boundary = compute_thicknesses(
        mat_for_thickness, mask=mask, spacing=spacing, layer_split=layer_split)
    base_th = base_thickness if base_thickness is not None else auto_base
    spin_th = spin_thickness if spin_thickness is not None else auto_spin

    spin_mask = mask[boundary:, :, :] if boundary is not None else mask
    base_mask = mask[:boundary, :, :] if boundary is not None else mask

    spin_mask_2d = spin_mask.any(axis=0)
    base_mask_2d = base_mask.any(axis=0)

    nx = xdim + 1
    ny = ydim + 1

    def node_id(ix, iy):
        return iy * nx + ix + 1

    nodes = []
    for iy in range(ny):
        y = oy + iy * sy
        for ix in range(nx):
            x = ox + ix * sx
            nodes.append((node_id(ix, iy), x, y, oz))

    elements = []
    base_eids = []
    spin_eids = []
    eid = 1
    for iy in range(ydim):
        for ix in range(xdim):
            if not base_mask_2d[iy, ix] and not spin_mask_2d[iy, ix]:
                continue
            n1 = node_id(ix, iy)
            n2 = node_id(ix + 1, iy)
            n3 = node_id(ix + 1, iy + 1)
            n4 = node_id(ix, iy + 1)
            elements.append((eid, (n1, n2, n3, n4)))
            if spin_mask_2d[iy, ix]:
                spin_eids.append(eid)
            else:
                base_eids.append(eid)
            eid += 1

    with open(path, 'w') as f:
        f.write('*Heading\n')
        f.write('** mat_to_shell_inp generated shell mesh\n')
        f.write('*Node\n')
        for nid, x, y, z in nodes:
            f.write('{:d}, {:.8e}, {:.8e}, {:.8e}\n'.format(nid, x, y, z))
        f.write('*Element, type=S4R\n')
        for eid, conn in elements:
            f.write('{:d}, {}\n'.format(eid, ', '.join(str(n) for n in conn)))
        if base_eids:
            _write_id_lines(f, base_eids, '*Elset, elset=BASE_SHELL')
        if spin_eids:
            _write_id_lines(f, spin_eids, '*Elset, elset=SPINODAL_SHELL')
        _write_id_lines(f, (nid for nid, _x, _y, _z in nodes), '*Nset, nset=ALLNODES')
        f.write('*Material, name={}\n'.format(material))
        if density is not None:
            f.write('*Density\n{:.6e}\n'.format(density))
        if elastic is not None:
            f.write('*Elastic\n{:.6e}, {:.6f}\n'.format(elastic[0], elastic[1]))
        f.write('*Shell Section, elset=BASE_SHELL, material={}, thickness={:.8e}\n'.format(
            material, base_th))
        if spin_eids:
            f.write('*Shell Section, elset=SPINODAL_SHELL, material={}, thickness={:.8e}\n'.format(
                material, base_th + spin_th))
        f.write('** End of file\n')

    return len(base_eids), len(spin_eids)


def convert_mat_to_shell(mat_path: str,
                         varname: Optional[str] = None,
                         spacing: Optional[float] = None,
                         origin: Optional[Tuple[float, float, float]] = None,
                         material: str = 'SPINODAL',
                         density: Optional[float] = None,
                         elastic: Optional[Tuple[float, float]] = None,
                         output_path: Optional[str] = None,
                         peel: int = 0,
                         base_thickness: Optional[float] = None,
                         spin_thickness: Optional[float] = None):
    mask, resolved_var, mat_data = load_mask_from_mat(mat_path, varname=varname, peel=peel)
    try:
        spacing_val = _resolve_spacing(spacing, mat_data)
    except ValueError:
        if 'voxelSizeXY' in mat_data:
            spacing_val = float(np.array(mat_data['voxelSizeXY']).astype(float).ravel()[0])
        else:
            raise
    origin_vec = _resolve_origin(origin, mat_data)
    layer_split = _detect_sheet_layers(mat_data)

    if output_path:
        out_path = output_path
    else:
        base, _ = os.path.splitext(mat_path)
        out_path = base + '_shell.inp'

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    base_count, spin_count = write_shell_inp(
        out_path,
        mask,
        spacing_val,
        origin_vec,
        material,
        density=density,
        elastic=elastic,
        mat_data=mat_data,
        layer_split=layer_split,
        base_thickness=base_thickness,
        spin_thickness=spin_thickness,
    )
    return out_path, base_count, spin_count


def main():
    args = parse_args()
    elastic_vals = tuple(args.elastic) if args.elastic is not None else None
    out_path, base_count, spin_count = convert_mat_to_shell(
        mat_path=args.mat,
        varname=args.var,
        spacing=args.spacing,
        origin=args.origin,
        material=args.material,
        density=args.density,
        elastic=elastic_vals,
        output_path=args.output,
        peel=args.peel,
        base_thickness=args.base_thickness,
        spin_thickness=args.spin_thickness,
    )
    print('Wrote shell mesh to {} (BASE elements: {}, SPIN elements: {})'.format(
        out_path, base_count, spin_count))


if __name__ == '__main__':
    main()
