%% Demo script for graph visualization
% Run this to visualize structural graphs

% Find paths automatically
scriptDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(scriptDir);  % .../Matlab
gnnRoot = fullfile(parentDir, 'GNN');

samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');
targetsDir = fullfile(gnnRoot, 'data', 'dataset', 'targets');

% Check if paths exist
if ~isfolder(samplesDir)
    error('Samples directory not found: %s', samplesDir);
end
if ~isfolder(targetsDir)
    error('Targets directory not found: %s', targetsDir);
end

fprintf('Samples Dir: %s\n', samplesDir);
fprintf('Targets Dir: %s\n\n', targetsDir);

%% Example 1: Single graph analysis

sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';  % Change sample ID here
options = struct('figName', 'single_graph_analysis', 'compareMode', 'none', 'showDeformation', false);

fprintf('Visualizing: %s\n', sample_id);
visualize_graph(sample_id, samplesDir, targetsDir, options);

fprintf('\n');

%% Example 2: Compare varying angles (comment/uncomment to try)

% options = struct('compareMode', 'angle', 'fixedValue', 1.0, ...
%                  'nSamples', 4, 'figName', 'varying_angle_comparison');
% visualize_graph('dummy', samplesDir, targetsDir, options);

%% Example 3: Compare varying thickness ratios (comment/uncomment to try)

% options = struct('compareMode', 'ratio', 'fixedValue', 30, ...
%                  'nSamples', 4, 'figName', 'varying_ratio_comparison');
% visualize_graph('dummy', samplesDir, targetsDir, options);

%% Example 4: Radius and neighbor ring topology

sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';  % Change sample ID here
options = struct('showNeighborRing', true, 'neighborNode', []);  % [] = auto-pick centroid node

fprintf('Visualizing radius/neighbor topology: %s\n', sample_id);
visualize_graph_topology(sample_id, samplesDir, options);