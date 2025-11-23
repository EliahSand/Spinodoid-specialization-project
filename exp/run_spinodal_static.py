#! /usr/bin/env python
"""Abaqus/CAE batch script to apply uniaxial loading to a voxel mesh."""


"""How to run:

abaqus cae noGUI=exp/run_spinodal_static.py -- -- Matlab/.../sheet.inp


"""


from abaqus import mdb, session
from abaqusConstants import *
import visualization
import odbAccess
import os
import sys


def parse_input_path():
    args = sys.argv
    if '--' in args:
        idx = args.index('--')
        cli_args = args[idx + 1:]
    else:
        cli_args = args[1:]
    if not cli_args:
        raise ValueError('Usage: abaqus cae noGUI=run_spinodal_static.py -- <mesh.inp>')
    inp_path = os.path.abspath(cli_args[0])
    if not os.path.isfile(inp_path):
        raise IOError('Input file not found: %s' % inp_path)
    return inp_path


def import_mesh(inp_path, model_name):
    if model_name in mdb.models:
        del mdb.models[model_name]
    model = mdb.ModelFromInputFile(name=model_name, inputFileName=inp_path)
    part_name = model.parts.keys()[0]
    part = model.parts[part_name]

    if 'SoftPolymer' not in model.materials:
        mat = model.Material(name='SoftPolymer')
        mat.Elastic(table=((1.0e6, 0.4),))
    if 'SpinodalSection' not in model.sections:
        model.HomogeneousSolidSection(name='SpinodalSection', material='SoftPolymer', thickness=None)
    region = part.Set(name='ALL_ELEMENTS', elements=part.elements)
    part.SectionAssignment(region=region, sectionName='SpinodalSection')

    assembly = model.rootAssembly
    assembly.DatumCsysByDefault(CARTESIAN)
    inst_name = assembly.instances.keys()[0]
    return model, assembly, inst_name


def create_step_and_bcs(model, assembly, inst_name):
    if 'LoadStep' not in model.steps:
        model.StaticStep(name='LoadStep', previous='Initial')

    # Ensure nodal outputs (U, S) are requested (explicitly nodal position)
    if 'NodalOutputs' in model.fieldOutputRequests:
        del model.fieldOutputRequests['NodalOutputs']
    model.FieldOutputRequest(name='NodalOutputs', createStepName='LoadStep',
                              variables=('U', 'S'))

    inst = assembly.instances[inst_name]
    xs = [n.coordinates[0] for n in inst.nodes]
    min_x = min(xs)
    max_x = max(xs)
    span = max(max_x - min_x, 1.0)
    tol = 1e-6 * span

    left_nodes = inst.nodes.getByBoundingBox(xMin=min_x - tol, xMax=min_x + tol,
                                             yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9)
    right_nodes = inst.nodes.getByBoundingBox(xMin=max_x - tol, xMax=max_x + tol,
                                              yMin=-1e9, yMax=1e9, zMin=-1e9, zMax=1e9)
    if not left_nodes:
        raise ValueError('No nodes detected on minimum-X face.')
    if not right_nodes:
        raise ValueError('No nodes detected on maximum-X face.')

    assembly.Set(name='BC_LEFT', nodes=left_nodes)
    assembly.Set(name='BC_RIGHT', nodes=right_nodes)

    if 'BC_Left' in model.boundaryConditions:
        del model.boundaryConditions['BC_Left']
    if 'BC_Right' in model.boundaryConditions:
        del model.boundaryConditions['BC_Right']

    model.EncastreBC(name='BC_Left', createStepName='Initial',
                     region=assembly.sets['BC_LEFT'])
    model.DisplacementBC(name='BC_Right', createStepName='LoadStep',
                         region=assembly.sets['BC_RIGHT'],
                         u1=0.1, u2=0.0, u3=0.0)


def run_job(model_name, job_name, work_dir):
    if job_name in mdb.jobs:
        del mdb.jobs[job_name]
    job = mdb.Job(name=job_name, model=model_name,
                  description='Spinodal tensile test',
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


def collect_exclusion_nodes(odb_root, instance_name):
    exclude = set()
    # Collect any sets whose name suggests spinodal layer
    spin_names = [n for n in odb_root.nodeSets.keys() if 'SPIN' in n.upper()]
    spin_names += [n for n in odb_root.elementSets.keys() if 'SPIN' in n.upper()]
    spin_names = list(set(spin_names))
    for name in spin_names:
        # Node sets
        if name in odb_root.nodeSets:
            try:
                nodes = odb_root.nodeSets[name].nodes[instance_name]
                for n in nodes:
                    exclude.add((instance_name, n.label))
            except KeyError:
                pass
        # Element sets
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
    odb = visualization.openOdb(path=odb_path)

    instance = odb.rootAssembly.instances.values()[0]
    inst_name = list(odb.rootAssembly.instances.keys())[0]
    xs = [node.coordinates[0] for node in instance.nodes]
    min_x = min(xs)
    max_x = max(xs)
    x_mid = 0.5 * (min_x + max_x)
    span = max(max_x - min_x, 1.0)
    unique_x = sorted(set(xs))
    dx_list = [b - a for a, b in zip(unique_x[:-1], unique_x[1:]) if b > a]
    if dx_list:
        min_dx = min(dx_list)
    else:
        min_dx = span
    slice_half = max(0.01 * min_dx, 1e-9)

    try:
        step = odb.steps['LoadStep']
    except KeyError:
        odb.close()
        raise ValueError('LoadStep not found in ODB.')
    frame = step.frames[-1]

    spinodal_nodes = collect_exclusion_nodes(odb.rootAssembly, inst_name)
    base_nodes = collect_base_nodes(odb.rootAssembly, inst_name)
    if base_nodes is not None:
        base_nodes = {n for n in base_nodes if n not in spinodal_nodes}

    # Compute a z-threshold from spinodal sets (exclude everything at/above min spinodal Z)
    spin_z_min = None
    spin_names = [n for n in odb.rootAssembly.nodeSets.keys() if 'SPIN' in n.upper()]
    spin_names += [n for n in odb.rootAssembly.elementSets.keys() if 'SPIN' in n.upper()]
    if spin_names:
        spin_node_labels = set()
        for name in spin_names:
            if name in odb.rootAssembly.nodeSets:
                try:
                    for n in odb.rootAssembly.nodeSets[name].nodes[inst_name]:
                        spin_node_labels.add(n.label)
                except KeyError:
                    pass
            if name in odb.rootAssembly.elementSets:
                for elem in odb.rootAssembly.elementSets[name].elements:
                    inst = getattr(elem, 'instanceName', inst_name)
                    if inst != inst_name:
                        continue
                    if hasattr(elem, 'connectivity'):
                        for lbl in elem.connectivity:
                            spin_node_labels.add(lbl)
        if spin_node_labels:
            spin_z = [instance.nodes[lbl-1].coordinates[2] for lbl in spin_node_labels if lbl <= len(instance.nodes)]
            if spin_z:
                spin_z_min = min(spin_z)

    u_field = frame.fieldOutputs['U']
    s_field = frame.fieldOutputs['S']

    # Map nodal displacements directly
    u_dict = {(getattr(v, 'instanceName', inst_name), v.nodeLabel): v.data for v in u_field.values}

    # Average element stresses to nodes
    elem_conn = {elem.label: list(elem.connectivity) for elem in instance.elements}
    stress_acc = {}
    stress_cnt = {}
    for val in s_field.values:
        conn = elem_conn.get(val.elementLabel, [])
        for nl in conn:
            key = (inst_name, nl)
            if key not in stress_acc:
                stress_acc[key] = [0.0]*len(val.data)
                stress_cnt[key] = 0
            stress_acc[key] = [a + b for a, b in zip(stress_acc[key], val.data)]
            stress_cnt[key] += 1

    s_dict = {}
    mises_dict = {}
    for key, acc in stress_acc.items():
        cnt = max(1, stress_cnt.get(key, 1))
        avg = [a/float(cnt) for a in acc]
        s_dict[key] = avg
        if len(avg) >= 6:
            s11, s22, s33, s12, s13, s23 = avg[:6]
            mises = ((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 + 6*(s12**2 + s13**2 + s23**2))**0.5 / (3**0.5)
        else:
            mises = 0.0
        mises_dict[key] = mises

    # Determine base z-levels; keep only the bottom surface of the base layer
    if base_nodes:
        base_z = [instance.nodes[node_label-1].coordinates[2] for (_inst, node_label) in base_nodes if _inst == inst_name]
        base_z_min = min(base_z) if base_z else None
    else:
        base_z_min = min(node.coordinates[2] for node in instance.nodes)

    rows = []
    for node in instance.nodes:
        key = (inst_name, node.label)
        if key in spinodal_nodes:
            continue
        if base_nodes is not None and key not in base_nodes:
            continue
        if spin_z_min is not None and node.coordinates[2] >= spin_z_min - 1e-9:
            continue
        if base_z_min is not None and node.coordinates[2] > base_z_min + 1e-9:
            continue
        x = node.coordinates[0]
        if abs(x - x_mid) > slice_half:
            continue
        disp = u_dict.get(key)
        stress = s_dict.get(key)
        mises = mises_dict.get(key, 0.0)
        if disp is None or stress is None:
            continue
        rows.append((node.label, node.coordinates, disp, stress, mises))

    output_path = os.path.join(fea_dir, 'midplane_results.csv')
    with open(output_path, 'w') as fh:
        fh.write('Label,X,Y,Z,U1,U2,U3,S11,S22,S33,S12,S13,S23,SMises\n')
        for label, coords, disp, stress, mises in rows:
            fh.write('{:d},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n'.format(
                label,
                coords[0], coords[1], coords[2],
                disp[0], disp[1], disp[2],
                stress[0], stress[1], stress[2],
                stress[3], stress[4], stress[5],
                mises))
    # Optional: write a minimal field report for displacements (safe position)
    rpt_path = os.path.join(fea_dir, 'midplane_results.rpt')
    try:
        session.writeFieldReport(fileName=rpt_path, append=OFF, sortItem='Node Label', odb=odb,
                                  step=0, frame=frame.incrementNumber,
                                  outputPosition=INTEGRATION_POINT, variable=(('U', INTEGRATION_POINT), ('S', INTEGRATION_POINT)))
    except Exception:
        # Ignore report errors so CSV generation still succeeds
        pass
    odb.close()
    return output_path


def main():
    inp_path = parse_input_path()
    work_dir = os.path.dirname(inp_path)
    fea_dir = os.path.join(work_dir, 'FEA')
    if not os.path.isdir(fea_dir):
        os.makedirs(fea_dir)

    model_name = 'SpinodalModel'
    model, assembly, inst_name = import_mesh(inp_path, model_name)
    create_step_and_bcs(model, assembly, inst_name)

    job_name = os.path.splitext(os.path.basename(inp_path))[0] + '_job'
    odb_path = run_job(model_name, job_name, fea_dir)
    csv_path = extract_midplane_results(odb_path, fea_dir)

    print('Analysis complete. Results stored in %s' % fea_dir)
    print('Mid-plane data: %s' % csv_path)


if __name__ == '__main__':
    main()
