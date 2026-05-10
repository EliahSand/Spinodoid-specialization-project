function [figModel, figFull] = visualize_graph(sample_id, samplesDir, varargin)
%VISUALIZE_GRAPH Show the two graph views that matter for one GNN sample.
%   Tab 1: reduced input graph exactly as the active model loader gives it
%   to the model: node features [x, y, radius, boundary] and edge_index for
%   hybrid graphs, [x, y, radius] for legacy graph-only data, or the
%   schema v4 MAT/spline control graph.
%
%   Tab 2: full reference graph for the same sample, rebuilt from the raw
%   sheet_shell.inp with the gnn_prep_spinodal fullGraph builder.
%
%   Name/value options:
%     ShowFullGraph    true by default. Set false for standalone sample.mat files.
%     ShowMatEnvelopeTab true by default for schema v4 graphs with debug masks.
%     ShowSplineCurves true by default for schema v4 graphs.
%     ShowRadiusDisks  true by default.
%
%   Extra legacy arguments are accepted and ignored.

if nargin < 2 || isempty(samplesDir)
    samplesDir = default_samples_dir();
end

sample_id = char(string(sample_id));
samplesDir = char(string(samplesDir));
opts = parse_visualize_options(varargin{:});

[gnnRoot, repoRoot] = infer_roots(samplesDir);
helpersDir = fullfile(gnnRoot, 'helpers');
if isfolder(helpersDir)
    addpath(helpersDir);
end
gnnGraphDir = fullfile(gnnRoot, 'datasetCreation', 'gnn_graph');
if isfolder(gnnGraphDir)
    addpath(gnnGraphDir);
end

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
if opts.ShowFullGraph && ~isfolder(gnnPrepRoot)
    error('visualize_graph:MissingGnnPrep', ...
        'Could not find gnn_prep_spinodal at: %s', gnnPrepRoot);
end
if opts.ShowFullGraph
    addpath(gnnPrepRoot);
    addpath(genpath(fullfile(gnnPrepRoot, 'src')));
end

[X, EI, featureNames, sampleData] = load_model_input_graph(sample_id, samplesDir);

assert_model_graph_schema(X, EI, sample_id);

figModel = figure('Name', sprintf('Graph views: %s', sample_id), ...
    'Color', 'w', 'Position', [80, 100, 900, 760]);
figFull = [];
tabs = uitabgroup(figModel);
tabModel = uitab(tabs, 'Title', 'GNN input');

axModel = axes('Parent', tabModel);
plot_model_input_graph(axModel, X, EI, featureNames, sampleData, opts);
title(axModel, sprintf('GNN input graph: %s', escape_title(sample_id)));

if opts.ShowMatEnvelopeTab && can_plot_mat_envelope(sampleData)
    tabEnvelope = uitab(tabs, 'Title', 'MAT envelope');
    axEnvelope = axes('Parent', tabEnvelope);
    plot_mat_envelope(axEnvelope, sampleData, opts);
end

if opts.ShowFullGraph
    figFull = figModel;
    tabFull = uitab(tabs, 'Title', 'Full graph');

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
end

sampleFigDir = fullfile(samplesDir, sample_id);
if isfolder(sampleFigDir)
    figPath = fullfile(sampleFigDir, sprintf('%s_graph_views.fig', sample_id));
    savefig(figModel, figPath);
    fprintf('Saved figure to %s\n', figPath);
end

fprintf('Sample: %s\n', sample_id);
fprintf('GNN input: %d nodes, %d edges, features [%s]\n', ...
    size(X, 2), size(EI, 2), feature_name_text(featureNames));
if isfield(sampleData, 'gnn_data') && isfield(sampleData.gnn_data, 'schema_version')
    fprintf('Schema: v%d', sampleData.gnn_data.schema_version);
    if isfield(sampleData.gnn_data, 'representation')
        fprintf(' (%s)', sampleData.gnn_data.representation);
    end
    if isfield(sampleData.gnn_data, 'branch_count')
        fprintf(', branches: %d', sampleData.gnn_data.branch_count);
    end
    fprintf('\n');
end
if opts.ShowFullGraph
    fprintf('Full graph: %d nodes, %d edges, source %s\n', ...
        fullGraph.num_nodes, size(fullGraph.edges_local, 1), inpPath);
end
end

function samplesDir = default_samples_dir()
scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset_hybrid', 'samples');
end

function [gnnRoot, repoRoot] = infer_roots(samplesDir)
datasetDir = fileparts(samplesDir);
dataDir = fileparts(datasetDir);
gnnRoot = fileparts(dataDir);
repoRoot = fileparts(fileparts(gnnRoot));
end

function [X, EI, featureNames, sampleData] = load_model_input_graph(sample_id, samplesDir)
sampleData = struct();
sampleMatPath = fullfile(samplesDir, sample_id, 'sample.mat');
if isfile(sampleMatPath)
    sampleData = load(sampleMatPath);
    sampleData.sample_mat_path = sampleMatPath;
    gd = sampleData.gnn_data;
    X = orient_feature_matrix(double(gd.x), get_field_default(gd, 'num_nodes', NaN));
    EI = double(gd.edge_index);
    featureNames = infer_feature_names(size(X, 1), gd);
    return;
end

if is_hybrid_samples_dir(samplesDir) || is_mat_spline_samples_dir(samplesDir)
    [X_cell, ei_cell] = load_hybrid_graph_dataset({sample_id}, samplesDir);
    X = double(X_cell{1});
    EI = double(ei_cell{1});
    featureNames = infer_feature_names(size(X, 1), struct());
else
    [X_cell, ei_cell] = load_graph_dataset({sample_id}, samplesDir);
    X = double(X_cell{1});
    EI = double(ei_cell{1});
    featureNames = {'x', 'y', 'radius'};
end
end

function tf = is_hybrid_samples_dir(samplesDir)
normalized = strrep(char(samplesDir), '\', filesep);
parts = strsplit(normalized, filesep);
tf = any(strcmp(parts, 'dataset_hybrid'));
end

function tf = is_mat_spline_samples_dir(samplesDir)
normalized = strrep(char(samplesDir), '\', filesep);
parts = strsplit(normalized, filesep);
tf = any(strcmp(parts, 'dataset_mat_spline'));
end

function assert_model_graph_schema(X, EI, sample_id)
if size(X, 1) < 3
    error('visualize_graph:BadModelFeatures', ...
        'Expected model features for "%s" to have at least 3 rows. Got %dx%d.', ...
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

function plot_model_input_graph(ax, X, EI, featureNames, sampleData, opts)
x = X(1, :).';
y = X(2, :).';
r = X(3, :).';
schemaVersion = get_schema_version(sampleData);
isMatSpline = schemaVersion == 4 || has_feature(featureNames, 'is_internal_control');

hold(ax, 'on');
draw_edges(ax, x, y, EI, [0.72, 0.72, 0.72], 0.45);

validDisks = isfinite(r) & r > 0;
if opts.ShowRadiusDisks
    if isMatSpline
        draw_radius_disks(ax, x(validDisks), y(validDisks), r(validDisks), ...
            [0.20, 0.42, 0.85], 0.16);
    else
        draw_radius_disks(ax, x(validDisks), y(validDisks), r(validDisks), ...
            [0.20, 0.42, 0.85], 0.92);
    end
end

if isMatSpline
    if opts.ShowSplineCurves
        draw_branch_splines(ax, X, sampleData);
    end
    endpoint = feature_logical(X, featureNames, 'is_endpoint', 5);
    junction = feature_logical(X, featureNames, 'is_junction', 6);
    internal = feature_logical(X, featureNames, 'is_internal_control', 7) & ~endpoint & ~junction;
    unmarked = ~(endpoint | junction | internal);

    if any(internal)
        scatter(ax, x(internal), y(internal), 18, [0.20, 0.42, 0.85], ...
            'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.45, ...
            'DisplayName', 'internal control');
    end
    if any(endpoint)
        scatter(ax, x(endpoint), y(endpoint), 30, [0.90, 0.20, 0.18], ...
            'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.55, ...
            'DisplayName', 'endpoint');
    end
    if any(junction)
        scatter(ax, x(junction), y(junction), 42, [0.05, 0.05, 0.05], ...
            's', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.55, ...
            'DisplayName', 'junction');
    end
    if any(unmarked)
        scatter(ax, x(unmarked), y(unmarked), 12, [0.20, 0.42, 0.85], ...
            'filled', 'MarkerEdgeColor', 'none');
    end
else
    boundary = feature_logical(X, featureNames, 'boundary', 4);
    if any(~validDisks)
        scatter(ax, x(~validDisks), y(~validDisks), 12, [0.20, 0.42, 0.85], ...
            'filled', 'MarkerEdgeColor', 'none');
    end
    if any(boundary)
        scatter(ax, x(boundary), y(boundary), 18, [0.90, 0.20, 0.18], ...
            'o', 'LineWidth', 0.8);
    end
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

function draw_radius_disks(ax, x, y, r, faceColor, faceAlpha)
if isempty(x)
    return;
end
if nargin < 5 || isempty(faceColor)
    faceColor = [0.20, 0.42, 0.85];
end
if nargin < 6 || isempty(faceAlpha)
    faceAlpha = 0.92;
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
    'FaceColor', faceColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', faceAlpha);
end

function draw_edges(ax, x, y, EI, color, lineWidth)
if isempty(EI)
    return;
end

edgeX = [x(EI(1, :)).'; x(EI(2, :)).'; nan(1, size(EI, 2))];
edgeY = [y(EI(1, :)).'; y(EI(2, :)).'; nan(1, size(EI, 2))];
plot(ax, edgeX(:), edgeY(:), '-', 'Color', color, 'LineWidth', lineWidth);
end

function draw_branch_splines(ax, X, sampleData)
if ~isfield(sampleData, 'matSplineGraph') || ~isfield(sampleData.matSplineGraph, 'branches')
    return;
end

branches = sampleData.matSplineGraph.branches;
if isempty(branches)
    return;
end

xy = X(1:2, :).';
curveColor = [0.05, 0.19, 0.39];
for b = 1:numel(branches)
    ctrlIds = double(branches(b).control_indices(:));
    ctrlIds = ctrlIds(ctrlIds >= 1 & ctrlIds <= size(xy, 1));
    if numel(ctrlIds) < 2
        continue;
    end
    ctrl = xy(ctrlIds, :);
    curve = evaluate_branch_curve(ctrl, branches(b));
    if isempty(curve)
        if isfield(branches(b), 'is_cycle') && branches(b).is_cycle && size(ctrl, 1) > 2
            ctrl = [ctrl; ctrl(1, :)]; %#ok<AGROW>
        end
        plot(ax, ctrl(:, 1), ctrl(:, 2), '-', 'Color', curveColor, 'LineWidth', 1.2);
    else
        plot(ax, curve(:, 1), curve(:, 2), '-', 'Color', curveColor, 'LineWidth', 1.4);
        if isfield(branches(b), 'is_cycle') && branches(b).is_cycle && ...
                norm(curve(end, :) - curve(1, :)) > 10 * eps
            plot(ax, [curve(end, 1), curve(1, 1)], [curve(end, 2), curve(1, 2)], ...
                '-', 'Color', curveColor, 'LineWidth', 1.4);
        end
    end
end
end

function curve = evaluate_branch_curve(ctrl, branch)
curve = [];
if ~isfield(branch, 'degree') || ~isfield(branch, 'knots')
    return;
end

nCtrl = size(ctrl, 1);
degree = double(branch.degree);
knots = double(branch.knots(:).');
if degree < 1 || nCtrl <= degree || numel(knots) ~= nCtrl + degree + 1
    return;
end

t = linspace(0, 1, max(40, 12 * nCtrl)).';
B = bspline_basis_matrix(t, nCtrl, degree, knots);
curve = B * ctrl;
end

function B = bspline_basis_matrix(t, nCtrl, degree, knots)
t = t(:);
B = zeros(numel(t), nCtrl);
for i = 1:nCtrl
    B(:, i) = bspline_basis(i, degree, t, knots, nCtrl);
end
B(t == 1, :) = 0;
B(t == 1, nCtrl) = 1;
end

function N = bspline_basis(i, degree, t, knots, nCtrl)
if degree == 0
    N = double(knots(i) <= t & t < knots(i + 1));
    if i == nCtrl
        N(t == 1) = 1;
    end
    return;
end

denomA = knots(i + degree) - knots(i);
denomB = knots(i + degree + 1) - knots(i + 1);
A = zeros(size(t));
B = zeros(size(t));
if denomA > eps
    A = ((t - knots(i)) / denomA) .* bspline_basis(i, degree - 1, t, knots, nCtrl);
end
if denomB > eps
    B = ((knots(i + degree + 1) - t) / denomB) .* bspline_basis(i + 1, degree - 1, t, knots, nCtrl);
end
N = A + B;
end

function tf = can_plot_mat_envelope(sampleData)
tf = isfield(sampleData, 'matSplineGraph') && ...
    isfield(sampleData, 'matDebug') && ...
    isfield(sampleData.matDebug, 'occupancy_mask') && ...
    isfield(sampleData.matDebug, 'x_grid') && ...
    isfield(sampleData.matDebug, 'y_grid');
end

function plot_mat_envelope(ax, sampleData, opts)
occMask = logical(sampleData.matDebug.occupancy_mask);
xVals = double(sampleData.matDebug.x_grid(:));
yVals = double(sampleData.matDebug.y_grid(:));

[reconMask, sampled, metrics] = reconstruct_mat_spline_graph( ...
    sampleData.matSplineGraph, xVals, yVals, ...
    'OccupancyMask', occMask, ...
    'SamplesPerGridStep', opts.EnvelopeSamplesPerGridStep);
fprintf('MAT envelope: sampled=%d IoU=%.4f missing=%.2f%% extra=%.2f%% one-sided=%.6g\n', ...
    sampled.sample_count, metrics.iou, 100 * metrics.missing_fraction, ...
    100 * metrics.extra_fraction, metrics.one_sided_missing_error);

rgb = reconstruction_rgb(occMask, reconMask);
image(ax, xVals, yVals, rgb);
set(ax, 'YDir', 'normal');
hold(ax, 'on');

draw_sampled_medial_circles(ax, sampled.points, opts.EnvelopeCircleStride);
draw_sampled_centerlines(ax, sampled.points, sampled.branch_index);

axis(ax, 'equal');
axis(ax, 'tight');
grid(ax, 'off');
box(ax, 'on');
xlabel(ax, 'x');
ylabel(ax, 'y');
title(ax, sprintf('MAT envelope: IoU %.3f, missing %.2f%%, extra %.2f%%, one-sided %.3g', ...
    metrics.iou, 100 * metrics.missing_fraction, 100 * metrics.extra_fraction, ...
    metrics.one_sided_missing_error));
set(ax, 'FontSize', 14, 'LineWidth', 1.0);
hold(ax, 'off');
end

function rgb = reconstruction_rgb(occMask, reconMask)
rgb = ones([size(occMask), 3]);

extra = reconMask & ~occMask;
covered = occMask & reconMask;
missing = occMask & ~reconMask;

rgb = paint_mask(rgb, extra, [1.00, 0.80, 0.38]);
rgb = paint_mask(rgb, covered, [0.69, 0.86, 1.00]);
rgb = paint_mask(rgb, missing, [0.92, 0.20, 0.18]);
end

function rgb = paint_mask(rgb, mask, color)
for c = 1:3
    channel = rgb(:, :, c);
    channel(mask) = color(c);
    rgb(:, :, c) = channel;
end
end

function draw_sampled_medial_circles(ax, points, stride)
if isempty(points)
    return;
end
stride = max(1, stride);
maxCircles = 1000;
idx = 1:stride:size(points, 1);
if numel(idx) > maxCircles
    idx = idx(round(linspace(1, numel(idx), maxCircles)));
end
points = points(idx, :);

theta = linspace(0, 2 * pi, 28);
for i = 1:size(points, 1)
    r = points(i, 3);
    if ~isfinite(r) || r <= 0
        continue;
    end
    cx = points(i, 1) + r * cos(theta);
    cy = points(i, 2) + r * sin(theta);
    plot(ax, cx, cy, '-', 'Color', [0.50, 0.78, 0.96], 'LineWidth', 0.35);
end
end

function draw_sampled_centerlines(ax, points, branchIndex)
if isempty(points)
    return;
end
ids = unique(branchIndex(:).');
for b = ids
    mask = branchIndex == b;
    if nnz(mask) < 2
        continue;
    end
    p = points(mask, :);
    plot(ax, p(:, 1), p(:, 2), '-', 'Color', [0.02, 0.12, 0.28], 'LineWidth', 1.0);
end
end

function opts = parse_visualize_options(varargin)
opts = struct();
opts.ShowFullGraph = true;
opts.ShowMatEnvelopeTab = true;
opts.ShowSplineCurves = true;
opts.ShowRadiusDisks = true;
opts.EnvelopeSamplesPerGridStep = 2.5;
opts.EnvelopeCircleStride = 3;

k = 1;
while k <= numel(varargin)
    key = varargin{k};
    if (ischar(key) || isstring(key)) && k < numel(varargin)
        value = varargin{k + 1};
        switch lower(char(string(key)))
            case 'showfullgraph'
                opts.ShowFullGraph = option_to_logical(value);
            case 'showmatenvelopetab'
                opts.ShowMatEnvelopeTab = option_to_logical(value);
            case 'showsplinecurves'
                opts.ShowSplineCurves = option_to_logical(value);
            case 'showradiusdisks'
                opts.ShowRadiusDisks = option_to_logical(value);
            case 'envelopesamplespergridstep'
                opts.EnvelopeSamplesPerGridStep = double(value);
            case 'envelopecirclestride'
                opts.EnvelopeCircleStride = max(1, round(double(value)));
        end
        k = k + 2;
    else
        k = k + 1;
    end
end
end

function tf = option_to_logical(value)
if islogical(value) || isnumeric(value)
    tf = logical(value);
elseif ischar(value) || isstring(value)
    tf = any(strcmpi(char(string(value)), {'true', 'on', 'yes', '1'}));
else
    tf = logical(value);
end
tf = tf(1);
end

function X = orient_feature_matrix(rawX, numNodes)
if ~isnan(numNodes)
    if size(rawX, 1) == numNodes
        X = rawX.';
        return;
    elseif size(rawX, 2) == numNodes
        X = rawX;
        return;
    end
end

if size(rawX, 1) <= 16 && size(rawX, 2) > size(rawX, 1)
    X = rawX;
else
    X = rawX.';
end
end

function value = get_field_default(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function featureNames = infer_feature_names(nFeatures, gd)
if isstruct(gd) && isfield(gd, 'feature_names') && ~isempty(gd.feature_names)
    featureNames = cellstr(string(gd.feature_names));
    return;
end

if isstruct(gd) && isfield(gd, 'schema_version') && gd.schema_version == 4
    featureNames = {'x', 'y', 'radius', 'degree', ...
        'is_endpoint', 'is_junction', 'is_internal_control'};
elseif nFeatures == 4
    featureNames = {'x', 'y', 'radius', 'boundary'};
elseif nFeatures == 3
    featureNames = {'x', 'y', 'radius'};
elseif nFeatures == 7
    featureNames = {'x', 'y', 'radius', 'degree', ...
        'is_endpoint', 'is_junction', 'is_internal_control'};
else
    featureNames = arrayfun(@(i) sprintf('f%d', i), 1:nFeatures, ...
        'UniformOutput', false);
end
end

function txt = feature_name_text(featureNames)
if isempty(featureNames)
    txt = '';
else
    txt = strjoin(cellstr(string(featureNames)), ' ');
end
end

function version = get_schema_version(sampleData)
version = NaN;
if isfield(sampleData, 'gnn_data') && isfield(sampleData.gnn_data, 'schema_version')
    version = double(sampleData.gnn_data.schema_version);
end
end

function tf = has_feature(featureNames, name)
tf = feature_row(featureNames, name, 0) > 0;
end

function values = feature_logical(X, featureNames, name, fallbackRow)
row = feature_row(featureNames, name, fallbackRow);
if row > 0 && size(X, 1) >= row
    values = logical(X(row, :).');
else
    values = false(size(X, 2), 1);
end
end

function row = feature_row(featureNames, name, fallbackRow)
row = 0;
if ~isempty(featureNames)
    names = lower(string(featureNames));
    hit = find(names == lower(string(name)), 1);
    if ~isempty(hit)
        row = hit;
        return;
    end
end
if nargin >= 3
    row = fallbackRow;
end
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
