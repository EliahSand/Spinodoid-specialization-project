% step9_aggregate_mat_spline_graphs.m
%
% Aggregates schema v4 MAT/spline control graphs into a single file.
% This writes under dataset_mat_spline and leaves dataset_hybrid v3 intact.

%% Config

if ~exist('REPO_ROOT', 'var')
    REPO_ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
if ~exist('SAMPLES_DIR', 'var')
    SAMPLES_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_mat_spline', 'samples');
end
if ~exist('TARGETS_DIR', 'var')
    TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_mat_spline', 'targets');
end

SCHEMA_VERSION = 4;
RASTER_SIZE = 128;

%% Discover sample folders

items = dir(SAMPLES_DIR);
sampleDirs = {};
for k = 1:numel(items)
    if items(k).isdir && items(k).name(1) ~= '.'
        matPath = fullfile(SAMPLES_DIR, items(k).name, 'sample.mat');
        if isfile(matPath)
            sampleDirs{end + 1, 1} = items(k).name; %#ok<AGROW>
        end
    end
end

n = numel(sampleDirs);
if n == 0
    error('step9_mat_spline:no_graphs', 'No sample.mat files found under %s', SAMPLES_DIR);
end
fprintf('Found %d MAT/spline v4 graphs. Aggregating...\n', n);

%% Load and aggregate

sample_ids = sampleDirs;
X_cell = cell(n, 1);
ei_cell = cell(n, 1);
N_vec = zeros(n, 1);
TR_vec = nan(n, 1);
ANG_vec = nan(n, 1);
Dense_cell = cell(n, 1);
has_dense_count = 0;
feature_names = {};
representation = 'mat_spline_control_graph';

t0 = tic;
for i = 1:n
    d = load(fullfile(SAMPLES_DIR, sample_ids{i}, 'sample.mat'));
    gd = d.gnn_data;
    if ~isfield(gd, 'schema_version') || gd.schema_version ~= SCHEMA_VERSION
        error('step9_mat_spline:bad_schema', ...
            'Expected schema_version=%d in %s, got %s.', ...
            SCHEMA_VERSION, sample_ids{i}, mat2str(get_field_default(gd, 'schema_version', NaN)));
    end

    X_cell{i} = single(gd.x.');
    ei_cell{i} = double(gd.edge_index);
    N_vec(i) = gd.num_nodes;
    TR_vec(i) = gd.tr_ratio;
    ANG_vec(i) = gd.ang_deg;
    if isempty(feature_names) && isfield(gd, 'feature_names')
        feature_names = gd.feature_names;
    end

    if isfield(d, 'dense_data') && isfield(d.dense_data, 'raster')
        Dense_cell{i} = uint8(d.dense_data.raster);
        has_dense_count = has_dense_count + 1;
    else
        Dense_cell{i} = rasterize_graph_xy(X_cell{i}, RASTER_SIZE);
    end

    if mod(i, 500) == 0
        fprintf('  %d / %d  (%.1f s)\n', i, n, toc(t0));
    end
end
fprintf('Loaded %d MAT/spline graphs in %.1f s\n', n, toc(t0));

%% Save aggregate

if ~isfolder(TARGETS_DIR), mkdir(TARGETS_DIR); end

created_at = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
schema_version = SCHEMA_VERSION;

outFile = fullfile(TARGETS_DIR, 'graphs_all.mat');
save(outFile, 'sample_ids', 'X_cell', 'ei_cell', 'N_vec', 'TR_vec', 'ANG_vec', ...
    'Dense_cell', 'feature_names', 'representation', 'schema_version', 'created_at', '-v7.3');
fprintf('Saved: %s\n', outFile);

%% Manifest

try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

manifest = struct( ...
    'n_graphs', n, ...
    'n_dense_native', has_dense_count, ...
    'raster_size', RASTER_SIZE, ...
    'schema_version', schema_version, ...
    'representation', representation, ...
    'feature_names', {feature_names}, ...
    'git_sha', git_sha, ...
    'created_at', created_at);

manifestFile = fullfile(TARGETS_DIR, 'graphs_all_manifest.json');
fid = fopen(manifestFile, 'w');
fprintf(fid, '%s\n', jsonencode(manifest));
fclose(fid);
fprintf('Manifest: %s\n', manifestFile);
fprintf('Done.\n');

function val = get_field_default(s, fieldName, defaultVal)
if isfield(s, fieldName)
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function raster = rasterize_graph_xy(X, gridSize)
xy = double(X(1:2, :).');
xLo = min(xy(:, 1)); xHi = max(xy(:, 1));
yLo = min(xy(:, 2)); yHi = max(xy(:, 2));
cols = 1 + round((xy(:, 1) - xLo) / max(xHi - xLo, eps) * (gridSize - 1));
rows = 1 + round((xy(:, 2) - yLo) / max(yHi - yLo, eps) * (gridSize - 1));
cols = min(gridSize, max(1, cols));
rows = min(gridSize, max(1, rows));
occupancy = false(gridSize, gridSize);
idx = sub2ind([gridSize, gridSize], rows, cols);
occupancy(idx) = true;
raster = uint8(cat(3, occupancy, false(gridSize, gridSize)));
end
