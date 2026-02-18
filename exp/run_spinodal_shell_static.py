#! /usr/bin/env python
"""Abaqus/CAE batch script to apply uniaxial loading on a spinodal shell mesh.

Sections:
  BaseShellSection     -> BASE_SHELL (base layer thickness)
  SpinodalShellSection -> SPINODAL_SHELL (base + spinodal thickness)
BCs:
  Left face encastre (min X); right face displacement u1=0.1*(maxX-minX).

Usage:
    abaqus cae noGUI=exp/run_spinodal_shell_static.py -- path/to/shell_mesh.inp 
"""

from abaqus import mdb, session
from abaqusConstants import *
import visualization  # noqa: F401
import odbAccess  # noqa: F401
import os
import sys
import json
import time

import numpy as np

try:
    import scipy.io as sio
except ImportError:  # pragma: no cover
    sio = None
try:  # pragma: no cover
    import h5py
except ImportError:  # pragma: no cover
    h5py = None


def parse_input_paths():
    args = sys.argv
    if '--' in args:
        idx = args.index('--')
        cli_args = args[idx + 1:]
    else:
        cli_args = args[1:]
    # Filter out known flags that Abaqus injects (e.g., "-cae"), but keep other tokens (they may belong to spaced paths).
    skip_flags = set(['-cae', '-mesa'])
    filtered = []
    i = 0
    while i < len(cli_args):
        p = cli_args[i]
        pl = p.lower()
        if pl in skip_flags:
            i += 1
            continue
        if pl == '-lmlog':
            i += 1  # skip the log path token if present
            if i < len(cli_args) and not cli_args[i].lower().endswith('.inp'):
                i += 1
            continue
        if pl.startswith('-nogui') or 'nogui=' in pl:
            i += 1
            continue
        if pl.startswith('-') and not os.path.isfile(p):
            i += 1
            continue
        filtered.append(p)
        i += 1
    cli_args = filtered
    if not cli_args:
        raise ValueError('Usage: abaqus cae noGUI=exp/run_spinodal_shell_static.py -- <mesh.inp> [...]')

    # Reconstruct paths that may have been split on spaces (e.g., OneDrive paths).
    inp_paths = []
    i = 0
    while i < len(cli_args):
        found = None
        for j in range(len(cli_args), i, -1):
            cand = ' '.join(cli_args[i:j])
            if cand.lower().endswith('.inp') and os.path.isfile(cand):
                found = cand
                i = j
                break
        if found:
            inp_paths.append(os.path.abspath(found))
            continue
        i += 1

    if not inp_paths:
        raise IOError('No valid input files found; checked: %s' % (' '.join(cli_args)))
    return inp_paths


def parse_manifest(inp_path):
    manifest_path = os.path.join(os.path.dirname(inp_path), 'mesh_manifest.json')
    if not os.path.exists(manifest_path):
        raise IOError('mesh_manifest.json not found next to INP')
    with open(manifest_path, 'r') as fh:
        manifest = json.load(fh)
    return manifest_path, manifest


def parse_runlog_thickness(log_path):
    """Extract t_base_mm and t_spin_mm from a run_log.txt if present."""
    if not os.path.isfile(log_path):
        return None, None
    base = None
    spin = None
    try:
        with open(log_path, 'r') as fh:
            for line in fh:
                line_low = line.strip().lower()
                if line_low.startswith('t_base_mm'):
                    try:
                        base = float(line.split(':', 1)[1].strip()) / 1000.0
                    except Exception:
                        pass
                elif line_low.startswith('t_spin_mm'):
                    try:
                        spin = float(line.split(':', 1)[1].strip()) / 1000.0
                    except Exception:
                        pass
    except Exception:
        return None, None
    return base, spin


def load_mat(mat_path):
    if sio is None:
        raise ImportError('scipy (or h5py for v7.3) required to read MAT files')
    try:
        data = sio.loadmat(mat_path)
    except NotImplementedError:
        if h5py is None:
            raise RuntimeError('MAT v7.3 detected. Install h5py or supply fallback %s' % mat_path)
        data = {}
        with h5py.File(mat_path, 'r') as f:
            for key in f.keys():
                data[key] = f[key][...]
    return data


def resolve_spacing(manifest, mat_data):
    if manifest and manifest.get('spacing') is not None:
        arr = np.array(manifest['spacing']).astype(float).ravel()
        if arr.size == 1:
            return float(arr[0])
        if arr.size == 3:
            return float(arr[0]), float(arr[1]), float(arr[2])
    for key in ('voxelSpacing', 'voxel_size', 'spacing', 'voxelSizeXY'):
        if key in mat_data:
            arr = np.array(mat_data[key]).astype(float).ravel()
            if arr.size == 1:
                return float(arr[0])
            if arr.size == 3:
                return float(arr[0]), float(arr[1]), float(arr[2])
    raise ValueError('Voxel spacing not found (manifest.spacing or voxelSpacing / voxelSizeXY missing).')


def resolve_origin(manifest, mat_data):
    if manifest and manifest.get('origin') is not None:
        arr = np.array(manifest['origin']).astype(float).ravel()
        if arr.size >= 3:
            return float(arr[0]), float(arr[1]), float(arr[2])
    if 'origin' in mat_data:
        arr = np.array(mat_data['origin']).astype(float).ravel()
        if arr.size >= 3:
            return float(arr[0]), float(arr[1]), float(arr[2])
    return (0.0, 0.0, 0.0)


def _normalize_spacing(spacing):
    arr = np.array(spacing, dtype=float).ravel()
    if arr.size == 1:
        return float(arr[0]), float(arr[0]), float(arr[0])
    if arr.size == 3:
        return float(arr[0]), float(arr[1]), float(arr[2])
    raise ValueError('Spacing must be scalar or length-3 sequence, got %s' % (arr,))


def resolve_mask(mat_data, var_name=None):
    chosen = None
    if var_name:
        if var_name not in mat_data:
            raise KeyError('Mask variable %s not in MAT (or any 3D mask)' % var_name)
        chosen = var_name
        arr = np.array(mat_data[var_name])
    else:
        arr = None
        for key, value in mat_data.items():
            if key.startswith('__'):
                continue
            cand = np.array(value)
            if cand.ndim == 3:
                arr = cand
                chosen = key
                break
    if arr is None:
        raise RuntimeError('No 3D mask found; please set "var" in manifest.')
    mask = np.squeeze(np.array(arr)).astype(bool)
    if mask.ndim != 3:
        raise ValueError('Mask must be 3D, got %s' % (mask.shape,))
    return mask, chosen


def infer_boundary(mat_data, mask):
    boundary = None
    if 'zVoxelThickness' in mat_data:
        zth = np.array(mat_data['zVoxelThickness']).astype(float).ravel()
        if zth.size >= 2:
            first = zth[0]
            tol = max(1e-12, 1e-6 * max(1.0, abs(first)))
            diff_idx = np.where(np.abs(zth - first) > tol)[0]
            if diff_idx.size > 0:
                boundary = int(diff_idx[0])
    if boundary is None:
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
    return max(0, min(int(boundary), mask.shape[0]))


def infer_thicknesses(mat_data, mask, spacing, boundary):
    _, _, sz = _normalize_spacing(spacing)
    base_th = 0.0
    spin_th = 0.0
    if 'zVoxelThickness' in mat_data:
        zth = np.array(mat_data['zVoxelThickness']).astype(float).ravel()
        if boundary <= zth.size:
            base_th = float(np.sum(zth[:boundary])) if boundary > 0 else 0.0
            spin_th = float(np.sum(zth[boundary:])) if boundary < zth.size else 0.0
    if base_th == 0.0 and spin_th == 0.0:
        base_th = boundary * sz
        spin_th = max(mask.shape[0] - boundary, 0) * sz
    return base_th, spin_th


def _ensure_material_with_elastic(model):
    """Return material name; create SoftPolymer with elastic if needed."""
    if 'SoftPolymer' in model.materials:
        mat_name = 'SoftPolymer'
    elif model.materials:
        mat_name = model.materials.keys()[0]
    else:
        mat_name = 'SoftPolymer'
        model.Material(name=mat_name)
    mat = model.materials[mat_name]
    if not hasattr(mat, 'elastic'):
        mat.Elastic(table=((1.0e6, 0.4),))
    return mat_name


def _deduce_thicknesses_from_sections(model):
    base_th = None
    spin_total = None
    fallback = None
    for name, sec in model.sections.items():
        if not hasattr(sec, 'thickness'):
            continue
        t = float(sec.thickness)
        if fallback is None:
            fallback = t
        if 'BASE' in name.upper():
            base_th = t
        if 'SPIN' in name.upper():
            spin_total = t
    return base_th, spin_total, fallback


def classify_elements_by_mask(part, mask2d, spacing, origin):
    sx, sy, _ = _normalize_spacing(spacing)
    ox, oy = origin[0], origin[1]
    ydim, xdim = mask2d.shape

    def to_index(coord, origin_val, step, dim):
        idx = int(round((coord - origin_val) / float(step) - 0.5))
        if idx < 0:
            return 0
        if idx >= dim:
            return dim - 1
        return idx

    base_eids = []
    spin_eids = []
    for elem in part.elements:
        conn = getattr(elem, 'connectivity', ())
        coords = [part.nodes[nid - 1].coordinates for nid in conn]
        if not coords:
            continue
        x = sum(p[0] for p in coords) / float(len(coords))
        y = sum(p[1] for p in coords) / float(len(coords))
        ix = to_index(x, ox, sx, xdim)
        iy = to_index(y, oy, sy, ydim)
        if mask2d[iy, ix]:
            spin_eids.append(elem.label)
        else:
            base_eids.append(elem.label)
    return base_eids, spin_eids


def _is_laminate_case(manifest, inp_path):
    text = ' '.join([
        str(manifest.get('notes', '')),
        str(manifest.get('mask', '')),
        str(manifest.get('mat', '')),
        str(manifest.get('var', '')),
        str(inp_path),
    ]).lower()
    return 'laminate' in text


def collect_exclusion_nodes(odb_root, instance_name):
    exclude = set()
    spin_names = [n for n in odb_root.nodeSets.keys() if 'SPIN' in n.upper()]
    spin_names += [n for n in odb_root.elementSets.keys() if 'SPIN' in n.upper()]
    spin_names = list(set(spin_names))
    for name in spin_names:
        if name in odb_root.nodeSets:
            try:
                nodes = odb_root.nodeSets[name].nodes[instance_name]
                for n in nodes:
                    exclude.add((instance_name, n.label))
            except KeyError:
                pass
        if name in odb_root.elementSets:
            elems = odb_root.elementSets[name].elements
            for elem in elems:
                inst = getattr(elem, 'instanceName', instance_name)
                if inst != instance_name:
                    continue
                if hasattr(elem, 'connectivity'):
                    for lbl in elem.connectivity:
                        exclude.add((inst, lbl))
    return exclude


def collect_base_nodes(odb_root, instance_name):
    include = set()
    base_names = [n for n in odb_root.elementSets.keys() if 'BASE' in n.upper()]
    if not base_names:
        return None
    for name in base_names:
        elems = odb_root.elementSets[name].elements
        for elem in elems:
            inst = getattr(elem, 'instanceName', instance_name)
            if inst != instance_name:
                continue
            if hasattr(elem, 'connectivity'):
                for lbl in elem.connectivity:
                    include.add((inst, lbl))
    return include if include else None


def extract_midplane_results(odb_path, fea_dir):
    """
    Post-process the ODB to mimic:
        'selection view of the middle nodes -> probe values -> write to file'
    in batch/noGUI mode.

    - Finds geometric mid-plane along X.
    - Selects nodes in a thin band around that plane.
    - Extracts U (nodal displacements) and S (stresses), computes von Mises.
    - Writes midplane_results_shell.csv in fea_dir.
    """
    rows = []
    odb = None
    try:
        from odbAccess import openOdb
        from abaqusConstants import NODAL

        odb = openOdb(path=odb_path)

        # --- Instance ---
        instances = list(odb.rootAssembly.instances.keys())
        if not instances:
            raise ValueError('No instances found in ODB.')
        inst_name = instances[0]
        instance = odb.rootAssembly.instances[inst_name]

        # --- Mid-plane in X ---
        xs = [node.coordinates[0] for node in instance.nodes]
        if not xs:
            raise ValueError('No nodes found in instance %s.' % inst_name)

        min_x = min(xs)
        max_x = max(xs)
        span = max(max_x - min_x, 1.0)

        unique_x = sorted(set(xs))
        dx_list = [b - a for a, b in zip(unique_x[:-1], unique_x[1:]) if b > a]
        min_dx = dx_list[0] if dx_list else span

        x_mid = 0.5 * (min_x + max_x)
        target_x = min(unique_x, key=lambda xv: abs(xv - x_mid)) if unique_x else x_mid
        slice_half = 0.5 * min_dx if min_dx > 0 else 0.5 * span

        # --- Frame (last frame, prefer 'LoadStep') ---
        frame = None
        if 'LoadStep' in odb.steps and odb.steps['LoadStep'].frames:
            frame = odb.steps['LoadStep'].frames[-1]
        else:
            for st in odb.steps.values():
                if st.frames:
                    frame = st.frames[-1]
                    break

        u_dict = {}
        s_dict = {}
        mises_dict = {}

        if frame is not None and 'U' in frame.fieldOutputs:
            u_field = frame.fieldOutputs['U']
            for v in u_field.values:
                key = (getattr(v, 'instanceName', inst_name), v.nodeLabel)
                u_dict[key] = tuple(v.data)

        if frame is not None and 'S' in frame.fieldOutputs:
            s_field = frame.fieldOutputs['S']

            # Try nodal stresses first
            try:
                s_nodal = s_field.getSubset(position=NODAL)
                s_values = list(s_nodal.values)
            except Exception:
                s_values = []

            if s_values:
                for val in s_values:
                    key = (getattr(val, 'instanceName', inst_name), val.nodeLabel)
                    data = list(val.data)
                    s_dict[key] = data

                    # Compute von Mises from whatever is available (pad to 6)
                    s = (data + [0.0] * 6)[:6]
                    s11, s22, s33, s12, s13, s23 = s
                    mises = ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2 +
                             6 * (s12 ** 2 + s13 ** 2 + s23 ** 2)) ** 0.5 / (3 ** 0.5)
                    mises_dict[key] = mises
            else:
                # Fallback: average IP stresses to nodes
                elem_conn = {elem.label: list(elem.connectivity) for elem in instance.elements}
                stress_acc = {}
                stress_cnt = {}

                for val in s_field.values:
                    if getattr(val, 'instanceName', inst_name) != inst_name:
                        continue
                    conn = elem_conn.get(val.elementLabel, [])
                    for nl in conn:
                        key = (inst_name, nl)
                        if key not in stress_acc:
                            stress_acc[key] = [0.0] * len(val.data)
                            stress_cnt[key] = 0
                        stress_acc[key] = [a + b for a, b in zip(stress_acc[key], val.data)]
                        stress_cnt[key] += 1

                for key, acc in stress_acc.items():
                    cnt = max(1, stress_cnt.get(key, 1))
                    avg = [a / float(cnt) for a in acc]
                    s_dict[key] = avg

                    s = (avg + [0.0] * 6)[:6]
                    s11, s22, s33, s12, s13, s23 = s
                    mises = ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2 +
                             6 * (s12 ** 2 + s13 ** 2 + s23 ** 2)) ** 0.5 / (3 ** 0.5)
                    mises_dict[key] = mises

        # --- Collect mid-plane nodes ---
        for node in instance.nodes:
            if slice_half > 0 and abs(node.coordinates[0] - target_x) > slice_half:
                continue

            key = (inst_name, node.label)
            disp = u_dict.get(key, (0.0, 0.0, 0.0))
            stress = s_dict.get(key, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            mises = mises_dict.get(key, 0.0)
            rows.append((node.label, node.coordinates, disp, stress, mises))

        # If nothing was captured, fall back to all nodes
        if not rows:
            for node in instance.nodes:
                key = (inst_name, node.label)
                disp = u_dict.get(key, (0.0, 0.0, 0.0))
                stress = s_dict.get(key, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                mises = mises_dict.get(key, 0.0)
                rows.append((node.label, node.coordinates, disp, stress, mises))

        print('Mid-plane extraction: total nodes=%d, selected=%d'
              % (len(instance.nodes), len(rows)))

    except Exception as exc:
        print('Post-processing warning (mid-plane extraction): %s' % exc)
    finally:
        try:
            if odb is not None:
                odb.close()
        except Exception:
            pass

    if not rows:
        # Last-resort dummy row so CSV is never empty
        rows.append((0, (0.0, 0.0, 0.0),
                     (0.0, 0.0, 0.0),
                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     0.0))

    # --- Write CSV safely ---
    if not os.path.isdir(fea_dir):
        os.makedirs(fea_dir)

    output_path = os.path.join(fea_dir, 'midplane_results_shell.csv')
    with open(output_path, 'w') as fh:
        fh.write('Label,X,Y,Z,U1,U2,U3,S_11,S_22,S_33,S_12,S_13,S_23,S_Mises\n')
        for label, coords, disp, stress, mises in rows:
            # Make everything safe-length
            x, y, z = (list(coords) + [0.0, 0.0, 0.0])[:3]
            u1, u2, u3 = (list(disp) + [0.0, 0.0, 0.0])[:3]
            s = (list(stress) + [0.0] * 6)[:6]
            s11, s22, s33, s12, s13, s23 = s

            fh.write(
                '{:d},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n'
                .format(
                    int(label),
                    float(x), float(y), float(z),
                    float(u1), float(u2), float(u3),
                    float(s11), float(s22), float(s33),
                    float(s12), float(s13), float(s23),
                    float(mises)
                )
            )

    print('Mid-plane probe-style data written to %s (rows=%d)'
          % (output_path, len(rows)))
    return output_path, len(rows)


def import_shell_mesh(inp_path, model_name, mask2d, spacing, origin, base_th, spin_th, preserve_input_sections=False):
    if model_name in mdb.models:
        del mdb.models[model_name]
    model = mdb.ModelFromInputFile(name=model_name, inputFileName=inp_path)
    part_name = model.parts.keys()[0]
    part = model.parts[part_name]
    if part.elements:
        elem_type = part.elements[0].type
        elem_type_str = str(elem_type)
        if not elem_type_str.upper().startswith('S'):
            raise ValueError('Input mesh must already be shell elements (S4/S3).')

    mat_name = _ensure_material_with_elastic(model)
    # Use supplied thicknesses if available; otherwise try to deduce from existing sections
    base_from_sec, spin_total_from_sec, fallback_th = _deduce_thicknesses_from_sections(model)
    if base_th is None:
        base_th = base_from_sec if base_from_sec is not None else (spin_total_from_sec if spin_total_from_sec is not None else fallback_th)
    if base_th is None:
        base_th = 1.0e-3
    if spin_th is None:
        if spin_total_from_sec is not None:
            spin_th = max(spin_total_from_sec - base_th, 1.0e-3)
        else:
            spin_th = 0.0

    base_section_th = base_th
    spin_section_th = base_th + spin_th

    # For laminate models, preserve section definitions/assignments from INP.
    if preserve_input_sections and len(part.sectionAssignments) > 0:
        print('Preserving shell sections from input INP (section assignments: %d).' % len(part.sectionAssignments))
        assembly = model.rootAssembly
        assembly.DatumCsysByDefault(CARTESIAN)
        inst_name = assembly.instances.keys()[0]
        return model, assembly, inst_name

    # Rebuild sections with guaranteed elastic material
    if 'BaseShellSection' in model.sections:
        del model.sections['BaseShellSection']
    if 'SpinodalShellSection' in model.sections:
        del model.sections['SpinodalShellSection']
    model.HomogeneousShellSection(
        name='BaseShellSection',
        material=mat_name,
        thickness=base_section_th,
        thicknessType=UNIFORM,
        preIntegrate=ON,
        poissonDefinition=DEFAULT,
        temperature=GRADIENT,
        integrationRule=SIMPSON,
        numIntPts=5,
    )
    model.HomogeneousShellSection(
        name='SpinodalShellSection',
        material=mat_name,
        thickness=spin_section_th,
        thicknessType=UNIFORM,
        preIntegrate=ON,
        poissonDefinition=DEFAULT,
        temperature=GRADIENT,
        integrationRule=SIMPSON,
        numIntPts=5,
    )

    # Build element lists
    base_eids = []
    spin_eids = []
    if mask2d is not None and spacing is not None and origin is not None:
        base_eids, spin_eids = classify_elements_by_mask(part, mask2d, spacing, origin)
        if not base_eids and not spin_eids:
            print('Mask-based classification failed; falling back to existing sets.')
    if not base_eids and not spin_eids:
        # Try existing sets in the part
        spin_sets = [s for s in part.sets.keys() if 'SPIN' in s.upper()]
        base_sets = [s for s in part.sets.keys() if 'BASE' in s.upper()]
        for sname in spin_sets:
            spin_eids.extend([e.label for e in part.sets[sname].elements])
        for sname in base_sets:
            base_eids.extend([e.label for e in part.sets[sname].elements])
    if not base_eids and not spin_eids:
        # Assign everything to base
        base_eids = [e.label for e in part.elements]

    if base_eids:
        part.SetFromElementLabels(name='BASE_SHELL', elementLabels=base_eids)
        part.SectionAssignment(
            region=part.sets['BASE_SHELL'],
            sectionName='BaseShellSection',
            offsetType=BOTTOM_SURFACE,
            thicknessAssignment=FROM_SECTION,
        )
    if spin_eids:
        part.SetFromElementLabels(name='SPINODAL_SHELL', elementLabels=spin_eids)
        part.SectionAssignment(
            region=part.sets['SPINODAL_SHELL'],
            sectionName='SpinodalShellSection',
            offsetType=BOTTOM_SURFACE,
            thicknessAssignment=FROM_SECTION,
        )

    # Ensure no element is left without a section
    assigned = set(base_eids + spin_eids)
    remaining = [e.label for e in part.elements if e.label not in assigned]
    if remaining:
        part.SetFromElementLabels(name='BASE_SHELL_REMAINDER', elementLabels=remaining)
        part.SectionAssignment(
            region=part.sets['BASE_SHELL_REMAINDER'],
            sectionName='BaseShellSection',
            offsetType=BOTTOM_SURFACE,
            thicknessAssignment=FROM_SECTION,
        )

    print('Section assignment: BASE=%d, SPIN=%d, REM=%d' % (len(base_eids), len(spin_eids), len(remaining)))

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    inst_name = assembly.instances.keys()[0]
    return model, assembly, inst_name


def create_step_and_bcs(model, assembly, inst_name):
    if 'LoadStep' not in model.steps:
        model.StaticStep(name='LoadStep', previous='Initial')

    if 'NodalOutputs' in model.fieldOutputRequests:
        del model.fieldOutputRequests['NodalOutputs']
    model.FieldOutputRequest(
        name='NodalOutputs',
        createStepName='LoadStep',
        variables=('U', 'S'),
    )

    inst = assembly.instances[inst_name]
    xs = [n.coordinates[0] for n in inst.nodes]
    min_x = min(xs)
    max_x = max(xs)
    L_raw = max_x - min_x
    if L_raw <= 0.0:
        raise ValueError('Invalid sheet length along X (max_x - min_x = %.6e).' % L_raw)
    L = L_raw
    u1_disp = 0.1 * L_raw

    # --- End "partition" bands defined by X-planes (robust for discrete meshes) ---
    # Band size as a fraction of the length (via discrete X-planes).
    # Override at runtime with environment variable BC_END_FRAC (e.g. 0.002 for 0.2%).
    end_frac = float(os.environ.get('BC_END_FRAC', '0.005'))  # default: 1%

    # Use unique X coordinates (planes/columns). This avoids oversized sets when spacing is coarse.
    unique_x = sorted(set(xs))
    nplanes = len(unique_x)
    if nplanes < 2:
        raise ValueError('Not enough unique X coordinates to define end bands (nplanes=%d).' % nplanes)

    # Number of X-planes to include on each end (at least 1)
    n_end = int(np.ceil(end_frac * nplanes))
    if n_end < 1:
        n_end = 1
    if n_end > nplanes:
        n_end = nplanes

    # Determine cut positions from the discrete planes
    x_left_max = unique_x[n_end - 1]
    x_right_min = unique_x[-n_end]

    # Small tolerance based on mesh extent
    tol = 1e-8 * max(1.0, abs(L))

    # Nodes in the left band: all nodes with x <= x_left_max
    left_nodes = inst.nodes.getByBoundingBox(
        xMin=min_x - tol, xMax=x_left_max + tol,
        yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9,
    )

    # Nodes in the right band: all nodes with x >= x_right_min
    right_nodes = inst.nodes.getByBoundingBox(
        xMin=x_right_min - tol, xMax=max_x + tol,
        yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9,
    )

    # Optional: element sets covering the same end bands
    left_elems = inst.elements.getByBoundingBox(
        xMin=min_x - tol, xMax=x_left_max + tol,
        yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9,
    )
    right_elems = inst.elements.getByBoundingBox(
        xMin=x_right_min - tol, xMax=max_x + tol,
        yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9,
    )

    print('BC band selection: end_frac=%.6g, nplanes=%d, n_end=%d, x_left_max=%.6e, x_right_min=%.6e, L=%.6e'
          % (end_frac, nplanes, n_end, x_left_max, x_right_min, L))

    if not left_nodes or not right_nodes:
        raise ValueError('Failed to find nodes on min/max X faces for BC sets.')

    assembly.Set(name='BC_LEFT', nodes=left_nodes)
    assembly.Set(name='BC_RIGHT', nodes=right_nodes)

    if left_elems:
        assembly.Set(name='BC_LEFT_ELEMS', elements=left_elems)
    if right_elems:
        assembly.Set(name='BC_RIGHT_ELEMS', elements=right_elems)

    if 'BC_Left' in model.boundaryConditions:
        del model.boundaryConditions['BC_Left']
    if 'BC_Right' in model.boundaryConditions:
        del model.boundaryConditions['BC_Right']

    model.EncastreBC(name='BC_Left', createStepName='Initial',
                     region=assembly.sets['BC_LEFT'])
    model.DisplacementBC(name='BC_Right', createStepName='LoadStep',
                         region=assembly.sets['BC_RIGHT'],
                         u1=u1_disp, u2=0.0, u3=0.0,
                         ur1=0.0,ur2 = 0.0, ur3=0.0
                         )


def run_job(model_name, job_name, work_dir):
    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    job = mdb.Job(name=job_name, model=model_name,
                  description='Spinodal shell tensile test',
                  numCpus=4, numDomains=4, multiprocessingMode=DEFAULT,
                  explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE)
    cwd = os.getcwd()
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    try:
        job.submit()
        job.waitForCompletion()
    finally:
        os.chdir(cwd)
    return os.path.join(work_dir, job_name + '.odb')


def main():
    inp_paths = parse_input_paths()
    for inp_path in inp_paths:
        manifest_path, manifest = parse_manifest(inp_path)
        mat_rel = manifest.get('mask') or manifest.get('mat')
        if not mat_rel:
            raise ValueError('Manifest missing "mask" or "mat" entry.')
        mat_path = os.path.join(os.path.dirname(manifest_path), mat_rel)
        mat_data = None
        mask = None
        mask2d = None
        spacing = None
        origin = None
        log_base, log_spin = parse_runlog_thickness(os.path.join(os.path.dirname(inp_path), 'run_log.txt'))
        base_th = log_base
        spin_th = log_spin
        if os.path.isfile(mat_path):
            try:
                mat_data = load_mat(mat_path)
                mask, mask_var = resolve_mask(mat_data, manifest.get('var'))
                spacing = resolve_spacing(manifest, mat_data)
                origin = resolve_origin(manifest, mat_data)
                boundary = infer_boundary(mat_data, mask)
                inf_base, inf_spin = infer_thicknesses(mat_data, mask, spacing, boundary)
                # Prefer run_log values when present
                if base_th is None:
                    base_th = inf_base
                if spin_th is None:
                    spin_th = inf_spin
                print('Using thicknesses: base = {:.6e} m, spinodal = {:.6e} m, spinodal section = {:.6e} m'.format(
                    base_th, spin_th, base_th + spin_th))
                spin_mask = mask[boundary:, :, :] if boundary is not None else mask
                mask2d = spin_mask.any(axis=0)
            except Exception as exc:
                print('MAT load failed (%s); falling back to existing INP sets/sections.' % exc)
        else:
            print('MAT file not found (%s); falling back to existing INP sets/sections.' % mat_path)

        work_dir = os.path.join(os.path.dirname(inp_path), 'FEA_shell')
        model_name = 'SpinodalShellModel'
        job_name = os.path.splitext(os.path.basename(inp_path))[0] + '_shell_job'

        preserve_input_sections = _is_laminate_case(manifest, inp_path)
        model, assembly, inst_name = import_shell_mesh(
            inp_path, model_name, mask2d, spacing, origin, base_th, spin_th,
            preserve_input_sections=preserve_input_sections)
        create_step_and_bcs(model, assembly, inst_name)
        # Save CAE snapshot alongside the shell FEA outputs.
        cae_path = os.path.join(work_dir, job_name + '.cae')
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir)
        mdb.saveAs(pathName=cae_path)
        odb_path = run_job(model_name, job_name, work_dir)
        print('Shell analysis complete. ODB: %s' % odb_path)
        try:
            csv_path, nrows = extract_midplane_results(odb_path, work_dir)
            print('Mid-plane data: %s (rows=%d)' % (csv_path, nrows))
        except Exception as exc:
            print('Post-processing skipped: %s' % exc)


if __name__ == '__main__':
    main()
