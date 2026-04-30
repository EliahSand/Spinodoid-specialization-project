# Spinodoid-Master
Author: Eliah Sand (NTNU, 2025)

End-to-end pipeline for predicting the midpoint deformation of a spinodoid sheet using hybrid CNN/GNN. 

## Requierements
- Matlab
- Abaqus
- Numpy
- Scipy
- hp5y



# Dataset Creation Pipeline

## Step 1 - Spinodoid generation

The spinodoid sheet is generated in MATLAB, with controllable parameters. This produces a **.mat** file

## Step 2 - mat to inp conversion

The .mat file is converted to a **.inp** file by a python script so it can be submitted for simulation in Abaqus

## Step 3 - graph generation

The spinodoid structure is captured and downsampled to a single **structural** graph, using a modified Zhang Zuen thinning algorithm. 

## Step 4 - Abaqus Simulation

The **.inp** file is sumbitted to Abaqus for simulation. Applied displacement on right side nodes, and fixed on the left end. *E=1e6* and *v=0.4* for all samples.



