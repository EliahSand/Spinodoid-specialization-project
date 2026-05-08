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


## Step 5 - Build u3 target dataset

Builds a target for the model for each sample as well as a aggregated **target.mat** file. 

## Step 6 - Split dataset

The dataset is split into training = 70, test = 15 and val = 15. This is done as a seperate step to ensure that the different model types are equally sampled across the different datasplits.

## Step 7 - Fit PCA on u3

PCA - Principal Component Analysis - is applied to the *training* dataset. It is only applied to this part to avoid data leakeage.

## step 8 - verify PCA target

This is just a precaution check, making sure that the PCA configurations are acceptable (to few PCA coefficients will not capture enough information)

## Step 9 - Aggreate graphs

Aggreates all strucutral graphs for the GNN as well as raster images for the CNN into a single file



