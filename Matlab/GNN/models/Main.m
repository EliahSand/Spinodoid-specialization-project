%% Graph neural network to predict the deformation of a spinodoid sheet

close all; clear; clc; 

rng(1); % Set random seed for reproducibility

%% Load dataset

%for training

load("Matlab/GNN/data/dataset/targets/pca_targets.mat", "Z_train");
load("Matlab/GNN/data/dataset/targets/split_indices.mat", "train_idx");

% for validation

load("GNN/data/dataset/targets/pca_targets.mat","Z_val");
load("GNN/data/dataset/targets/pca_model.mat","coeff","u3_mean");
load("GNN/data/dataset/targets/u3_targets.mat","U3_mat");
load("GNN/data/dataset/targets/split_indices.mat","val_idx");



%% 