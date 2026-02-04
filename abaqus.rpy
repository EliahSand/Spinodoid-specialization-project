# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022 replay file
# Internal Version: 2021_09_15-19.57.30 176069
# Run by eliahms on Wed Feb  4 14:17:48 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.99537, 1.02632), width=146.519, 
    height=101.811)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('exp/run_spinodal_shell_static.py', __main__.__dict__)
#: MAT load failed (MAT v7.3 detected. Install h5py or supply fallback C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\results\sheets\lamellar\sheetCone_tr50_ang045_lamellar_N64_1x20_run02\sheet.mat); falling back to existing INP sets/sections.
#: The model "SpinodalShellModel" has been created.
#: The part "PART-1" has been imported from the input file.
#: AbaqusException: Use of thickness modulus requires value for Poisson's ratio. This occurred while processing section keyword *SHELLSECTION. The keyword will be ignored. 
#: AbaqusException: Use of thickness modulus requires value for Poisson's ratio. This occurred while processing section keyword *SHELLSECTION. The keyword will be ignored. 
#: The model "SpinodalShellModel" has been imported from an input file. 
#: Please scroll up to check for error and warning messages.
#: Section assignment: BASE=40860, SPIN=41060, REM=0
#: The model database has been saved to "C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\results\sheets\lamellar\sheetCone_tr50_ang045_lamellar_N64_1x20_run02\FEA_shell\sheet_shell_shell_job.cae".
#: Shell analysis complete. ODB: C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\results\sheets\lamellar\sheetCone_tr50_ang045_lamellar_N64_1x20_run02\FEA_shell\sheet_shell_shell_job.odb
#: Model: C:/Users/eliahms/OneDrive - NTNU/Master 2026/Spinodoid-specialization-project/Matlab/results/sheets/lamellar/sheetCone_tr50_ang045_lamellar_N64_1x20_run02/FEA_shell/sheet_shell_shell_job.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       5
#: Number of Node Sets:          4
#: Number of Steps:              1
#: Mid-plane extraction: total nodes=83265, selected=65
#: Mid-plane probe-style data written to C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\results\sheets\lamellar\sheetCone_tr50_ang045_lamellar_N64_1x20_run02\FEA_shell\midplane_results_shell.csv (rows=65)
#: Mid-plane data: C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\results\sheets\lamellar\sheetCone_tr50_ang045_lamellar_N64_1x20_run02\FEA_shell\midplane_results_shell.csv (rows=65)
print 'RT script done'
#: RT script done
