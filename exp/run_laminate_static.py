#! /usr/bin/env python
"""Laminate solid Abaqus runner.

Usage:
    abaqus cae noGUI=exp/run_laminate_static.py -- -- path/to/sheet.inp
"""

import os
import sys

def _resolve_exp_dir():
    candidates = []
    if '__file__' in globals():
        candidates.append(os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) > 0 and sys.argv[0]:
        candidates.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    candidates.append(os.path.join(os.getcwd(), 'exp'))
    candidates.append(os.getcwd())
    for c in candidates:
        if c and os.path.isfile(os.path.join(c, 'run_spinodal_static.py')):
            return c
    return candidates[0]


THIS_DIR = _resolve_exp_dir()
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from run_spinodal_static import main


if __name__ == '__main__':
    main()
