%% Deprecated compatibility script
% The graph visualization no longer compares parameter sweeps. It shows one
% sample in the only two views needed for GNN debugging.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');

sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';
visualize_graph(sample_id, samplesDir);
