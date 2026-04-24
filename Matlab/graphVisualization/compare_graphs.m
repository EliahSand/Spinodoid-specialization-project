%% Compare graphs with different parameters

% Add graphVisualization to path
graphVizDir = fileparts(mfilename('fullpath'));
addpath(graphVizDir);

% Set paths (adjust to your setup)
samplesDir = fullfile(fileparts(fileparts(graphVizDir)), 'GNN', 'data', 'dataset', 'samples');
targetsDir = fullfile(fileparts(fileparts(graphVizDir)), 'GNN', 'data', 'dataset', 'targets');

% Choose comparison mode: 'angle' or 'ratio'
compareMode = 'angle';  % or 'ratio'

% Choose fixed value
% For angle mode: fix tr_ratio (e.g., 1.0 means 'tr100')
% For ratio mode: fix ang_deg (e.g., 30 means 'ang030')
fixedIndex = 100;  % tr=1.0 corresponds to 'tr100'

% Convert to string format
if strcmp(compareMode, 'angle')
    fixedValue = 'tr100';  % thickness ratio = 1.0
    options = struct('compareMode', 'angle', 'fixedValue', 1.0, ...
                     'nSamples', 4, 'figName', 'varying_angle_comparison');
else
    fixedValue = 'ang030';  % angle = 30 degrees
    options = struct('compareMode', 'ratio', 'fixedValue', 30, ...
                     'nSamples', 4, 'figName', 'varying_ratio_comparison');
end

% Use a dummy sample_id that contains the fixed parameter
% The script will filter automatically
sample_id = 'sheetCone_tr100_ang030_lamellar_N128_1x1';

% Visualize
visualize_graph(sample_id, samplesDir, targetsDir, options);