#!/usr/bin/env python3
"""Convert a cylindrical spinodoid stent mask (.mat) into an Abaqus .inp mesh.

This preserves the original cylindrical geometry by mapping the R/Theta/Z
voxel grid onto 3D coordinates before writing C3D8 elements.
"""

import argparse
import math
import os
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np

from mat_to_inp import (  # reuse MAT loading + origin resolver
    _resolve_origin,
    _write_id_lines,
    load_mask_from_mat,
)


def parse_args():
    p = argparse.ArgumentParser(
        description="MATLAB cylindrical stent mask -> Abaqus INP (curved geometry)",
    )
    p.add_argument("--mat", required=True, help="Path to .mat file (expects maskR/maskTheta/maskZ)")
    p.add_argument(
        "--var",
        default="solidMaskTiled",
        help="Variable name inside the .mat file (default: solidMaskTiled)",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output .inp path (default: <mat stem>_cyl.inp)",
    )
    p.add_argument("--material", default="SPINODAL", help="Material name for *Material card")
    p.add_argument("--density", type=float, default=None, help="Optional density value")
    p.add_argument(
        "--elastic",
        type=float,
        nargs=2,
        metavar=("E", "NU"),
        default=None,
        help="Optional elastic pair (E, nu)",
    )
    p.add_argument(
        "--origin",
        type=float,
        nargs=3,
        default=None,
        help="XYZ offset applied to node coordinates (default: from .mat or [0,0,0])",
    )
    p.add_argument(
        "--peel",
        type=int,
        default=0,
        help="Number of voxels to strip from each face before meshing (default 0)",
    )
    return p.parse_args()


def _as_1d_float(arr: np.ndarray, name: str) -> np.ndarray:
    vec = np.array(arr, dtype=float).ravel()
    if vec.size == 0:
        raise ValueError(f"{name} is empty")
    return vec


def _compute_edges(centers: np.ndarray, periodic: bool) -> np.ndarray:
    """Infer voxel boundary coordinates from center samples."""
    centers = np.array(centers, dtype=float).ravel()
    if centers.size == 0:
        raise ValueError("Center vector empty")
    if centers.size == 1:
        step = 1.0
    else:
        step = float(np.median(np.diff(centers)))

    if periodic:
        if centers.size < 2:
            raise ValueError("Need at least two theta samples for periodic axis")
        edges = np.empty_like(centers)
        edges[0] = centers[0] - 0.5 * step
        edges[1:] = 0.5 * (centers[:-1] + centers[1:])
        return edges

    edges = np.empty(centers.size + 1, dtype=float)
    edges[0] = centers[0] - 0.5 * step
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[-1] = centers[-1] + 0.5 * step
    return edges


def _orient_mask_rzt(mask: np.ndarray, r_len: int, z_len: int, t_len: int) -> np.ndarray:
    """Reorder mask axes to [r, z, theta] using length matches."""
    shape = mask.shape
    candidates = {
        "r": [i for i, d in enumerate(shape) if d == r_len],
        "z": [i for i, d in enumerate(shape) if d == z_len],
        "t": [i for i, d in enumerate(shape) if d == t_len],
    }
    if any(len(v) == 0 for v in candidates.values()):
        raise ValueError(
            f"Cannot align mask axes (shape {shape}) with lengths R={r_len}, Z={z_len}, Theta={t_len}"
        )
    # Prefer unique matches; otherwise pick the first unused axis per label.
    used = set()
    order = []
    for label in ("r", "z", "t"):
        axes = [ax for ax in candidates[label] if ax not in used]
        if not axes:
            axes = candidates[label]
        chosen = axes[0]
        used.add(chosen)
        order.append(chosen)
    return np.transpose(mask, axes=order)


def _node_id(r_idx: int, t_idx: int, z_idx: int, r_count: int, t_count: int) -> int:
    return z_idx * (t_count * r_count) + t_idx * r_count + r_idx + 1


def _write_nodes(
    fh,
    r_edges: Sequence[float],
    theta_edges: Sequence[float],
    z_edges: Sequence[float],
    origin: Tuple[float, float, float],
):
    r_edges = list(r_edges)
    theta_edges = list(theta_edges)
    z_edges = list(z_edges)
    ox, oy, oz = origin
    r_count = len(r_edges)
    t_count = len(theta_edges)

    node_idx = 1
    for z_val in z_edges:
        zc = oz + z_val
        for theta in theta_edges:
            ct = math.cos(theta)
            st = math.sin(theta)
            for r_val in r_edges:
                x = ox + r_val * ct
                y = oy + r_val * st
                fh.write(f"{node_idx:d}, {x:.8e}, {y:.8e}, {zc:.8e}\n")
                node_idx += 1
    return node_idx - 1  # total nodes


def _write_elements(
    fh,
    mask_rzt: np.ndarray,
    r_edges: Sequence[float],
    theta_edges: Sequence[float],
    z_edges: Sequence[float],
) -> Tuple[int, Iterable[int]]:
    r_cells, z_cells, t_cells = mask_rzt.shape
    r_count = len(r_edges)
    t_count = len(theta_edges)
    element_ids = []
    elem_id = 1

    for zi in range(z_cells):
        for ti in range(t_cells):
            t_next = (ti + 1) % t_cells  # wrap seam to close the tube
            for ri in range(r_cells):
                if not mask_rzt[ri, zi, ti]:
                    continue
                n1 = _node_id(ri, t_idx=ti, z_idx=zi, r_count=r_count, t_count=t_count)
                n2 = _node_id(ri + 1, t_idx=ti, z_idx=zi, r_count=r_count, t_count=t_count)
                n3 = _node_id(ri + 1, t_idx=t_next, z_idx=zi, r_count=r_count, t_count=t_count)
                n4 = _node_id(ri, t_idx=t_next, z_idx=zi, r_count=r_count, t_count=t_count)
                n5 = _node_id(ri, t_idx=ti, z_idx=zi + 1, r_count=r_count, t_count=t_count)
                n6 = _node_id(ri + 1, t_idx=ti, z_idx=zi + 1, r_count=r_count, t_count=t_count)
                n7 = _node_id(ri + 1, t_idx=t_next, z_idx=zi + 1, r_count=r_count, t_count=t_count)
                n8 = _node_id(ri, t_idx=t_next, z_idx=zi + 1, r_count=r_count, t_count=t_count)
                fh.write(f"{elem_id:d}, {n1}, {n2}, {n3}, {n4}, {n5}, {n6}, {n7}, {n8}\n")
                element_ids.append(elem_id)
                elem_id += 1
    return elem_id - 1, element_ids


def write_cyl_inp(
    path: str,
    mask_rzt: np.ndarray,
    r_centers: np.ndarray,
    z_centers: np.ndarray,
    theta_centers: np.ndarray,
    origin: Tuple[float, float, float],
    material: str,
    density: Optional[float],
    elastic: Optional[Tuple[float, float]],
):
    mask_rzt = mask_rzt.astype(bool)
    r_edges = _compute_edges(r_centers, periodic=False)
    r_edges[0] = max(0.0, r_edges[0])
    z_edges = _compute_edges(z_centers, periodic=False)
    theta_edges = _compute_edges(theta_centers, periodic=True)

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    with open(path, "w") as fh:
        fh.write("*Heading\n")
        fh.write("** mat_to_cyl_inp generated mesh (cylindrical)\n")
        fh.write("*Node\n")
        total_nodes = _write_nodes(fh, r_edges, theta_edges, z_edges, origin)

        fh.write("*Element, type=C3D8\n")
        total_elems, elem_ids = _write_elements(fh, mask_rzt, r_edges, theta_edges, z_edges)

        _write_id_lines(fh, elem_ids, "*Elset, elset=SOLID")
        _write_id_lines(fh, range(1, total_nodes + 1), "*Nset, nset=ALLNODES")

        fh.write(f"*Material, name={material}\n")
        if density is not None:
            fh.write(f"*Density\n{density:.6e}\n")
        if elastic is not None:
            fh.write(f"*Elastic\n{elastic[0]:.6e}, {elastic[1]:.6f}\n")
        fh.write(f"*Solid Section, elset=SOLID, material={material}\n,\n")
        fh.write("** End of file\n")

    return total_nodes, total_elems


def convert_cyl_mat_to_inp(
    mat_path: str,
    varname: Optional[str] = None,
    output_path: Optional[str] = None,
    material: str = "SPINODAL",
    density: Optional[float] = None,
    elastic: Optional[Tuple[float, float]] = None,
    origin: Optional[Tuple[float, float, float]] = None,
    peel: int = 0,
) -> Tuple[str, int]:
    mask, resolved_var, mat_data = load_mask_from_mat(mat_path, varname=varname, peel=peel)
    for key in ("maskR", "maskZ", "maskTheta"):
        if key not in mat_data:
            raise KeyError(f"{key} missing from MAT-file; required for cylindrical mapping.")
    r_centers = _as_1d_float(mat_data["maskR"], "maskR")
    z_centers = _as_1d_float(mat_data["maskZ"], "maskZ")
    theta_centers = _as_1d_float(mat_data["maskTheta"], "maskTheta")

    mask_rzt = _orient_mask_rzt(mask, len(r_centers), len(z_centers), len(theta_centers))
    origin_vec = _resolve_origin(origin, mat_data)

    if output_path:
        out_path = output_path
    else:
        base, ext = os.path.splitext(mat_path)
        out_path = f"{base}_cyl.inp"

    total_nodes, total_elems = write_cyl_inp(
        out_path,
        mask_rzt,
        r_centers,
        z_centers,
        theta_centers,
        origin_vec,
        material,
        density,
        elastic,
    )
    return out_path, int(mask_rzt.sum())


def main():
    args = parse_args()
    out_path, elem_count = convert_cyl_mat_to_inp(
        mat_path=args.mat,
        varname=args.var,
        output_path=args.output,
        material=args.material,
        density=args.density,
        elastic=tuple(args.elastic) if args.elastic is not None else None,
        origin=tuple(args.origin) if args.origin is not None else None,
        peel=args.peel,
    )
    print(f"Wrote cylindrical voxelized mesh with {elem_count} elements to {out_path}")


if __name__ == "__main__":
    main()
