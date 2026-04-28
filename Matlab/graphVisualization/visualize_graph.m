function [figModel, figFull] = visualize_graph(sample_id, samplesDir, varargin)
%VISUALIZE_GRAPH Show the two graph views that matter for one GNN sample.
%   Tab 1: reduced input graph exactly as load_graph_dataset gives it to the
%   model: node features [x, y, radius] and edge_index.
%
%   Tab 2: full reference graph for the same sample, rebuilt from the raw
%   sheet_shell.inp with the gnn_prep_spinodal fullGraph builder.
%
%   Extra legacy arguments are accepted and ignored.

if nargin < 2 || isempty(samplesDir)
    samplesDir = default_samples_dir();
end

sample_id = char(string(sample_id));
samplesDir = char(string(samplesDir));

[gnnRoot, repoRoot] = infer_roots(samplesDir);
addpath(fullfile(gnnRoot, 'helpers'));

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
if ~isfolder(gnnPrepRoot)
    error('visualize_graph:MissingGnnPrep', ...
        'Could not find gnn_prep_spinodal at: %s', gnnPrepRoot);
end
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

[X_cell, ei_cell] = load_graph_dataset({sample_id}, samplesDir);
X = double(X_cell{1});
EI = double(ei_cell{1});

assert_model_graph_schema(X, EI, sample_id);

figModel = figure('Name', sprintf('Graph views: %s', sample_id), ...
    'Color', 'w', 'Position', [80, 100, 900, 760]);
figFull = figModel;
tabs = uitabgroup(figModel);
tabModel = uitab(tabs, 'Title', 'GNN input');
tabFull = uitab(tabs, 'Title', 'Full graph');

axModel = axes('Parent', tabModel);
plot_model_input_graph(axModel, X, EI);
title(axModel, sprintf('GNN input graph: %s', escape_title(sample_id)));

inpPath = find_sample_inp(gnnRoot, sample_id);
inpData = read_abaqus_inp(inpPath);
fullGraph = build_full_reference_graph(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});

axFull = axes('Parent', tabFull);
plot_graph_mesh(fullGraph, ...
    'AxesHandle', axFull, ...
    'Title', sprintf('Full graph: %s', escape_title(sample_id)), ...
    'Style', 'paper2d', ...
    'BoundaryMode', 'interior', ...
    'XLabel', 'x', ...
    'YLabel', 'y');

fprintf('Sample: %s\n', sample_id);
fprintf('GNN input: %d nodes, %d edges, features [x y radius]\n', ...
    size(X, 2), size(EI, 2));
fprintf('Full graph: %d nodes, %d edges, source %s\n', ...
    fullGraph.num_nodes, size(fullGraph.edges_local, 1), inpPath);
end

function samplesDir = default_samples_dir()
scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');
end

function [gnnRoot, repoRoot] = infer_roots(samplesDir)
datasetDir = fileparts(samplesDir);
dataDir = fileparts(datasetDir);
gnnRoot = fileparts(dataDir);
repoRoot = fileparts(fileparts(gnnRoot));
end

function assert_model_graph_schema(X, EI, sample_id)
if size(X, 1) ~= 3
    error('visualize_graph:BadModelFeatures', ...
        'Expected model features for "%s" to be 3xN: [x; y; radius]. Got %dx%d.', ...
        sample_id, size(X, 1), size(X, 2));
end
if size(EI, 1) ~= 2
    error('visualize_graph:BadEdgeIndex', ...
        'Expected edge_index for "%s" to be 2xE. Got %dx%d.', ...
        sample_id, size(EI, 1), size(EI, 2));
end
if ~isempty(EI) && (min(EI(:)) < 1 || max(EI(:)) > size(X, 2))
    error('visualize_graph:EdgeIndexOutOfRange', ...
        'edge_index for "%s" references nodes outside 1..%d.', sample_id, size(X, 2));
end
end

function plot_model_input_graph(ax, X, EI)
x = X(1, :).';
y = X(2, :).';
r = X(3, :).';

hold(ax, 'on');
draw_edges(ax, x, y, EI, [0.72, 0.72, 0.72], 0.45);

validDisks = isfinite(r) & r > 0;
draw_radius_disks(ax, x(validDisks), y(validDisks), r(validDisks));
if any(~validDisks)
    scatter(ax, x(~validDisks), y(~validDisks), 12, [0.20, 0.42, 0.85], ...
        'filled', 'MarkerEdgeColor', 'none');
end

axis(ax, 'equal');
axis(ax, 'tight');
grid(ax, 'off');
box(ax, 'on');
xlabel(ax, 'x');
ylabel(ax, 'y');
set(ax, 'FontSize', 14, 'LineWidth', 1.0);
hold(ax, 'off');
end

function draw_radius_disks(ax, x, y, r)
if isempty(x)
    return;
end

nNodes = numel(x);
nSeg = 40;
theta = linspace(0, 2 * pi, nSeg + 1);
theta(end) = [];

vertices = zeros(nNodes * nSeg, 2);
faces = reshape(1:(nNodes * nSeg), nSeg, nNodes).';

for i = 1:nNodes
    rows = (i - 1) * nSeg + (1:nSeg);
    vertices(rows, 1) = x(i) + r(i) * cos(theta(:));
    vertices(rows, 2) = y(i) + r(i) * sin(theta(:));
end

patch(ax, ...
    'Faces', faces, ...
    'Vertices', vertices, ...
    'FaceColor', [0.20, 0.42, 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.92);
end

function draw_edges(ax, x, y, EI, color, lineWidth)
if isempty(EI)
    return;
end

edgeX = [x(EI(1, :)).'; x(EI(2, :)).'; nan(1, size(EI, 2))];
edgeY = [y(EI(1, :)).'; y(EI(2, :)).'; nan(1, size(EI, 2))];
plot(ax, edgeX(:), edgeY(:), '-', 'Color', color, 'LineWidth', lineWidth);
end

function inpPath = find_sample_inp(gnnRoot, sample_id)
trToken = regexp(sample_id, 'tr\d+', 'match', 'once');
angToken = regexp(sample_id, 'ang\d+', 'match', 'once');
rawSamplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');

if ~isempty(trToken) && ~isempty(angToken)
    inpPath = fullfile(rawSamplesRoot, trToken, angToken, sample_id, 'sheet_shell.inp');
    if isfile(inpPath)
        return;
    end
end

inpPath = find_inp_recursive(rawSamplesRoot, sample_id);
if isempty(inpPath)
    error('visualize_graph:MissingInp', ...
        'Could not find sheet_shell.inp for sample "%s" under %s.', ...
        sample_id, rawSamplesRoot);
end
end

function inpPath = find_inp_recursive(rootDir, sample_id)
inpPath = '';
if ~isfolder(rootDir)
    return;
end

items = dir(rootDir);
for k = 1:numel(items)
    name = items(k).name;
    if ~items(k).isdir || strcmp(name, '.') || strcmp(name, '..')
        continue;
    end

    subDir = fullfile(rootDir, name);
    if strcmp(name, sample_id)
        candidate = fullfile(subDir, 'sheet_shell.inp');
        if isfile(candidate)
            inpPath = candidate;
            return;
        end
    end

    inpPath = find_inp_recursive(subDir, sample_id);
    if ~isempty(inpPath)
        return;
    end
end
end

function out = escape_title(text)
out = strrep(char(string(text)), '_', '\_');
end
