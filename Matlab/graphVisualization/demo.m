%% Demo: visualize one GNN graph sample
% Produces exactly two figures:
%   1. the reduced input graph exactly as the GNN loader provides it
%   2. the matching full reference graph from sheet_shell.inp

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');

if ~isfolder(samplesDir)
    error('Samples directory not found: %s', samplesDir);
end

sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';
visualize_graph(sample_id, samplesDir);
