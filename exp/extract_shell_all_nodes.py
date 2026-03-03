#! /usr/bin/env python
"""Run shell analysis exactly like run_spinodal_shell_static, but export all nodes.

Usage:
    abaqus cae noGUI=exp/extract_shell_all_nodes.py -- path/to/shell_mesh.inp
"""

import os
import sys

from abaqusConstants import NODAL


def _resolve_exp_dir():
    candidates = []
    if '__file__' in globals():
        candidates.append(os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) > 0 and sys.argv[0]:
        candidates.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    candidates.append(os.path.join(os.getcwd(), 'exp'))
    candidates.append(os.getcwd())
    for c in candidates:
        if c and os.path.isfile(os.path.join(c, 'run_spinodal_shell_static.py')):
            return c
    return candidates[0]


THIS_DIR = _resolve_exp_dir()
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

import run_spinodal_shell_static as shell_runner  # noqa: E402


def extract_all_nodes_results(odb_path, fea_dir):
    """Drop-in replacement for shell_runner.extract_midplane_results."""
    rows = []
    odb = None
    try:
        # Reuse ready/open helpers when available.
        if hasattr(shell_runner, '_wait_for_odb_ready'):
            shell_runner._wait_for_odb_ready(odb_path)

        if hasattr(shell_runner, '_open_odb_with_retry'):
            opened = shell_runner._open_odb_with_retry(odb_path)
            if isinstance(opened, tuple):
                odb = opened[0]
            else:
                odb = opened
        else:
            from odbAccess import openOdb
            odb = openOdb(path=odb_path)

        instances = list(odb.rootAssembly.instances.keys())
        if not instances:
            raise ValueError('No instances found in ODB.')
        inst_name = instances[0]
        instance = odb.rootAssembly.instances[inst_name]

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
            try:
                s_nodal = s_field.getSubset(position=NODAL)
                s_values = list(s_nodal.values)
            except Exception:
                s_values = []

            if s_values:
                stress_acc = {}
                stress_cnt = {}
                for val in s_values:
                    key = (getattr(val, 'instanceName', inst_name), val.nodeLabel)
                    data = list(val.data)
                    if key not in stress_acc:
                        stress_acc[key] = [0.0] * len(data)
                        stress_cnt[key] = 0
                    stress_acc[key] = [a + b for a, b in zip(stress_acc[key], data)]
                    stress_cnt[key] += 1
                for key, acc in stress_acc.items():
                    cnt = max(1, stress_cnt.get(key, 1))
                    avg = [a / float(cnt) for a in acc]
                    s = (avg + [0.0] * 6)[:6]
                    s11, s22, s33, s12, s13, s23 = s
                    mises = ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2 +
                             6 * (s12 ** 2 + s13 ** 2 + s23 ** 2)) ** 0.5 / (3 ** 0.5)
                    s_dict[key] = s
                    mises_dict[key] = mises
            else:
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
                    s = (avg + [0.0] * 6)[:6]
                    s11, s22, s33, s12, s13, s23 = s
                    mises = ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2 +
                             6 * (s12 ** 2 + s13 ** 2 + s23 ** 2)) ** 0.5 / (3 ** 0.5)
                    s_dict[key] = s
                    mises_dict[key] = mises

        for node in instance.nodes:
            key = (inst_name, node.label)
            disp = u_dict.get(key, (0.0, 0.0, 0.0))
            stress = s_dict.get(key, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            mises = mises_dict.get(key, 0.0)
            rows.append((node.label, node.coordinates, disp, stress, mises))

        if not rows:
            rows.append((0, (0.0, 0.0, 0.0),
                        (0.0, 0.0, 0.0),
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        0.0))

        if not os.path.isdir(fea_dir):
            os.makedirs(fea_dir)
        output_path = os.path.join(fea_dir, 'midplane_results_shell.csv')
        with open(output_path, 'w') as fh:
            fh.write('Label,X,Y,Z,U1,U2,U3,S_11,S_22,S_33,S_12,S_13,S_23,S_Mises\n')
            for label, coords, disp, stress, mises in rows:
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

        print('All-node shell data written to %s (rows=%d)' % (output_path, len(rows)))
        return output_path, len(rows)
    finally:
        try:
            if odb is not None:
                odb.close()
        except Exception:
            pass


def main():
    # Keep full runner behavior identical; only swap extraction step.
    shell_runner.extract_midplane_results = extract_all_nodes_results
    shell_runner.main()


if __name__ == '__main__':
    main()
