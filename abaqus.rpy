# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022 replay file
# Internal Version: 2021_09_15-19.57.30 176069
# Run by eliahms on Tue May 12 13:23:06 2026
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
execfile('exp/run_spinodal_static.py', __main__.__dict__)
#: The model "SpinodalModel" has been created.
#: The part "PART-1" has been imported from the input file.
#: The model "SpinodalModel" has been imported from an input file. 
#: Please scroll up to check for error and warning messages.
#: The model database has been saved to "C:\Users\eliahms\OneDrive - NTNU\Master 2026\Spinodoid-specialization-project\Matlab\defectPrediction\results\lam045\crack_mid_edge_x\FEA\sheet_job.cae".
