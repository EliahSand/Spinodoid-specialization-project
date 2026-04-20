#!/usr/bin/env python
"""Run defect-only solid Abaqus jobs via the existing batch runner.

Usage:
    abaqus python Matlab/defectPrediction/run_batch_defect_abaqus.py [--overwrite] [--dry-run]

Equivalent direct command:
    abaqus python exp/batch_run_spinodal_solid.py --roots Matlab/defectPrediction/results
"""

from __future__ import print_function

import glob
import os
import sys


DEFECT_ROOT_REL = os.path.join('Matlab', 'defectPrediction', 'results')


def repo_root_from_here():
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.abspath(os.path.join(here, '..', '..'))


def format_extra_args(args):
    if not args:
        return ''
    return ' ' + ' '.join(args)


def main():
    repo_root = repo_root_from_here()
    defect_root = os.path.join(repo_root, DEFECT_ROOT_REL)
    manifests = []
    for dirpath, _dirnames, filenames in os.walk(defect_root):
        if 'mesh_manifest.json' in filenames:
            manifests.append(os.path.join(dirpath, 'mesh_manifest.json'))
    manifests.sort()
    if not manifests:
        raise SystemExit(
            'No generated defect manifests were found under %s.' % DEFECT_ROOT_REL
        )

    batch_script = os.path.join(repo_root, 'exp', 'batch_run_spinodal_solid.py')
    exp_dir = os.path.dirname(batch_script)
    if exp_dir not in sys.path:
        sys.path.insert(0, exp_dir)

    import batch_run_spinodal_solid as batch

    forwarded_args = sys.argv[1:]
    old_argv = list(sys.argv)
    old_cwd = os.getcwd()

    print('Running solid Abaqus batch for defect cases under %s' % DEFECT_ROOT_REL)
    print('[CMD ] abaqus python exp/batch_run_spinodal_solid.py --roots %s%s'
          % (DEFECT_ROOT_REL, format_extra_args(forwarded_args)))

    try:
        sys.argv = [batch_script, '--roots', DEFECT_ROOT_REL] + forwarded_args
        os.chdir(repo_root)
        batch.main()
    finally:
        sys.argv = old_argv
        os.chdir(old_cwd)


if __name__ == '__main__':
    main()
