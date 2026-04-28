% visualize_graph_topology.m
%
% Compatibility script for the simplified graph visualization.
% It shows exactly two figures:
%   1. the graph as load_graph_dataset gives it to the GNN
%   2. the matching full reference graph from sheet_shell.inp

scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fileparts(scriptDir);
matlabRoot = fileparts(gnnRoot);

addpath(fullfile(matlabRoot, 'graphVisualization'));

samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');
sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';

visualize_graph(sample_id, samplesDir);
