function [figModel, figFull] = visualize_mat_spline_graph(sample_id, samplesDir, varargin)
%VISUALIZE_MAT_SPLINE_GRAPH Visualize a schema v4 MAT/spline GNN sample.
%   Uses the same tabbed graphVisualization style as visualize_graph, but
%   defaults to Matlab/GNN/data/dataset_mat_spline/samples.

if nargin < 2 || isempty(samplesDir)
    samplesDir = default_mat_spline_samples_dir();
end

[figModel, figFull] = visualize_graph(sample_id, samplesDir, varargin{:});
end

function samplesDir = default_mat_spline_samples_dir()
scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset_mat_spline', 'samples');
end
