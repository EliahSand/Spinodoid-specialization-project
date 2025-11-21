"""Convert a MATLAB .mat file containing a 3D logical mask into an Abaqus .inp mesh.

Usage:
    python exp/mat_to_inp.py --mat solidMask.mat --var solidMask \
        --spacing 3.90625e-6 --peel 3 --output exp/results/spinodal_mesh.inp

Requirements: numpy + scipy (and h5py if the .mat file is v7.3 HDF5).
"""

import argparse
import os
from typing import Any, Dict, Optional, Tuple

import numpy as np

try:
    import scipy.io as sio
except ImportError:  # pragma: no cover
    raise ImportError('scipy is required to read MATLAB .mat files.')

try:  # pragma: no cover
    import h5py
except ImportError:  # pragma: no cover
    h5py = None


def parse_args():
    p = argparse.ArgumentParser(description='MATLAB mask -> Abaqus INP converter')
    p.add_argument('--mat', required=True, help='Path to .mat file containing the mask')
    p.add_argument('--var', default=None,
                   help='Variable name inside the .mat file (default: auto-detect first 3D logical array)')
    p.add_argument('--spacing', type=float, default=None,
                   help='Voxel edge length (meters). If omitted, attempts to read "voxelSpacing" from the .mat')
    p.add_argument('--output', default=None, help='Output .inp path (default: mat stem + .inp)')
    p.add_argument('--material', default='SPINODAL', help='Material name for *Material card')
    p.add_argument('--density', type=float, default=None, help='Optional density value')
    p.add_argument('--elastic', type=float, nargs=2, metavar=('E', 'NU'), default=None,
                   help='Optional elastic pair (E, nu)')
    p.add_argument('--origin', type=float, nargs=3, default=None,
                   help='XYZ offset applied to node coordinates (default: from .mat or [0,0,0])')
    p.add_argument('--peel', type=int, default=0,
                   help='Number of voxels to strip from each face before meshing (default 0)')
    return p.parse_args()


MAT_META_KEYS = {'__header__', '__version__', '__globals__'}


def _load_matfile(path: str) -> Dict[str, Any]:
    try:
        data = sio.loadmat(path)
    except NotImplementedError:
        if h5py is None:
            raise RuntimeError('MAT-file appears to be v7.3 (HDF5). Install h5py to read it.')
        data = {}
        with h5py.File(path, 'r') as f:
            for key in f.keys():
                data[key] = f[key][...]
    return data


def _choose_mask_variable(data: Dict[str, Any], requested: Optional[str]) -> Tuple[str, np.ndarray]:
    if requested:
        if requested not in data:
            raise KeyError('Variable %s not found' % requested)
        return requested, np.array(data[requested])
    candidates = []
    for key, value in data.items():
        if key in MAT_META_KEYS:
            continue
        arr = np.array(value)
        if arr.ndim == 3:
            candidates.append((key, arr))
    if not candidates:
        raise RuntimeError('No 3D variables found in MAT-file; please specify --var.')
    return candidates[0]


def load_mask_from_mat(path: str, varname: Optional[str] = None, peel: int = 0):
    data = _load_matfile(path)
    chosen_var, mask = _choose_mask_variable(data, varname)
    mask = np.array(mask)
    if mask.ndim == 0:
        raise ValueError('Mask variable %s is scalar' % chosen_var)
    mask = np.squeeze(mask).astype(bool)
    if mask.ndim != 3:
        raise ValueError('Mask must be 3D, got shape %s' % (mask.shape,))
    if not mask.any():
        raise ValueError('Mask contains no True voxels')
    mask = np.array(mask, order='C')
    if peel > 0:
        if peel*2 >= mask.shape[0] or peel*2 >= mask.shape[1] or peel*2 >= mask.shape[2]:
            raise ValueError('Peel value %d too large for mask shape %s' % (peel, mask.shape))
        mask[:peel, :, :] = False
        mask[-peel:, :, :] = False
        mask[:, :peel, :] = False
        mask[:, -peel:, :] = False
        mask[:, :, :peel] = False
        mask[:, :, -peel:] = False
    return mask, chosen_var, data


def _normalize_spacing(spacing):
    arr = np.array(spacing, dtype=float).ravel()
    if arr.size == 1:
        return float(arr[0]), float(arr[0]), float(arr[0])
    if arr.size == 3:
        return float(arr[0]), float(arr[1]), float(arr[2])
    raise ValueError('Spacing must be scalar or length-3 sequence, got %s' % (arr,))


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


def write_inp(path, mask, spacing, origin, material, density=None, elastic=None, layer_split=None):
    zdim, ydim, xdim = mask.shape
    ox, oy, oz = origin
    sx, sy, sz = _normalize_spacing(spacing)

    nx = xdim + 1
    ny = ydim + 1
    nz = zdim + 1

    def node_id(ix, iy, iz):
        return iz * (nx * ny) + iy * nx + ix + 1

    nodes = []
    for iz in range(nz):
        z = oz + iz * sz
        for iy in range(ny):
            y = oy + iy * sy
            for ix in range(nx):
                x = ox + ix * sx
                nodes.append((node_id(ix, iy, iz), x, y, z))

    elements = []
    layer_sets = None
    if layer_split and layer_split.get('axis', 0) == 0:
        boundary = int(layer_split.get('boundary', 0))
        layer_sets = {'BASE': [], 'SPIN': []}
    elem_id = 1
    for iz in range(zdim):
        for iy in range(ydim):
            for ix in range(xdim):
                if not mask[iz, iy, ix]:
                    continue
                n1 = node_id(ix,     iy,     iz)
                n2 = node_id(ix + 1, iy,     iz)
                n3 = node_id(ix + 1, iy + 1, iz)
                n4 = node_id(ix,     iy + 1, iz)
                n5 = node_id(ix,     iy,     iz + 1)
                n6 = node_id(ix + 1, iy,     iz + 1)
                n7 = node_id(ix + 1, iy + 1, iz + 1)
                n8 = node_id(ix,     iy + 1, iz + 1)
                elements.append((elem_id, (n1, n2, n3, n4, n5, n6, n7, n8)))
                if layer_sets is not None:
                    target = 'BASE' if iz < boundary else 'SPIN'
                    layer_sets[target].append(elem_id)
                elem_id += 1

    with open(path, 'w') as f:
        f.write('*Heading\n')
        f.write('** mat_to_inp generated mesh\n')
        f.write('*Node\n')
        for nid, x, y, z in nodes:
            f.write('{:d}, {:.8e}, {:.8e}, {:.8e}\n'.format(nid, x, y, z))
        f.write('*Element, type=C3D8\n')
        for eid, conn in elements:
            f.write('{:d}, {}\n'.format(eid, ', '.join(str(n) for n in conn)))
        _write_id_lines(f, (eid for eid, _ in elements), '*Elset, elset=SOLID')
        if layer_sets is not None:
            for name, ids in layer_sets.items():
                if not ids:
                    continue
                _write_id_lines(f, ids, '*Elset, elset={}'.format(name))
        _write_id_lines(f, (nid for nid, _x, _y, _z in nodes), '*Nset, nset=ALLNODES')
        f.write('*Material, name={}\n'.format(material))
        if density is not None:
            f.write('*Density\n{:.6e}\n'.format(density))
        if elastic is not None:
            f.write('*Elastic\n{:.6e}, {:.6f}\n'.format(elastic[0], elastic[1]))
        f.write('*Solid Section, elset=SOLID, material={}\n,\n'.format(material))
        f.write('** End of file\n')


def main():
    args = parse_args()
    out_path, elem_count = convert_mat_to_inp(
        mat_path=args.mat,
        varname=args.var,
        spacing=args.spacing,
        origin=args.origin,
        material=args.material,
        density=args.density,
        elastic=tuple(args.elastic) if args.elastic is not None else None,
        output_path=args.output,
        peel=args.peel,
    )
    print('Wrote voxelized mesh with {} elements to {}'.format(elem_count, out_path))


def _resolve_spacing(explicit: Optional[Any], mat_data: Dict[str, Any]):
    if explicit is not None:
        arr = np.array(explicit).astype(float).ravel()
        if arr.size == 1:
            return float(arr[0])
        if arr.size == 3:
            return float(arr[0]), float(arr[1]), float(arr[2])
        raise ValueError('Spacing override must be scalar or length-3.')
    for key in ('voxelSpacing', 'voxel_size', 'spacing'):
        if key in mat_data:
            val = np.array(mat_data[key]).astype(float).ravel()
            if val.size == 1:
                return float(val[0])
            if val.size == 3:
                return float(val[0]), float(val[1]), float(val[2])
    raise ValueError('Voxel spacing not provided; use --spacing or include "voxelSpacing" in the MAT-file.')


def _resolve_origin(explicit: Optional[Tuple[float, float, float]], mat_data: Dict[str, Any]) -> Tuple[float, float, float]:
    if explicit is not None:
        return tuple(float(v) for v in explicit)
    if 'origin' in mat_data:
        arr = np.array(mat_data['origin']).astype(float).ravel()
        if arr.size == 3:
            return float(arr[0]), float(arr[1]), float(arr[2])
    return (0.0, 0.0, 0.0)


def convert_mat_to_inp(mat_path: str,
                       varname: Optional[str] = None,
                       spacing: Optional[float] = None,
                       origin: Optional[Tuple[float, float, float]] = None,
                       material: str = 'SPINODAL',
                       density: Optional[float] = None,
                       elastic: Optional[Tuple[float, float]] = None,
                       output_path: Optional[str] = None,
                       peel: int = 0) -> Tuple[str, int]:
    mask, resolved_var, mat_data = load_mask_from_mat(mat_path, varname=varname, peel=peel)
    spacing_val = _resolve_spacing(spacing, mat_data)
    origin_vec = _resolve_origin(origin, mat_data)
    layer_split = _detect_sheet_layers(mat_data)
    if output_path:
        out_path = output_path
    else:
        base, _ = os.path.splitext(mat_path)
        out_path = base + '.inp'
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    write_inp(out_path, mask, spacing_val, origin_vec, material,
              density=density, elastic=elastic, layer_split=layer_split)
    return out_path, int(mask.sum())


def _detect_sheet_layers(mat_data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
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
    return {'axis': 0, 'boundary': boundary}


if __name__ == '__main__':
    main()
