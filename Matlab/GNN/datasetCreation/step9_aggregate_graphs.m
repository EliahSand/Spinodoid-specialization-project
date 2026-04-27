% step9_aggregate_graphs.m
%
% Aggregates all per-sample hybrid structural+raster graphs into a single
% file without touching the baseline data/dataset graph aggregate.
%
% Input:   data/dataset_hybrid/samples/<run>/sample.mat
% Outputs: data/dataset_hybrid/targets/graphs_all.mat
%          data/dataset_hybrid/targets/graphs_all_manifest.json
%
% Run after step 3 and before training.

%% Config

REPO_ROOT   = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SAMPLES_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_hybrid', 'samples');
TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_hybrid', 'targets');

SCHEMA_VERSION = 3;
RASTER_SIZE    = 128;

%% Discover sample folders

items = dir(SAMPLES_DIR);
sampleDirs = {};
for k = 1:numel(items)
    if items(k).isdir && items(k).name(1) ~= '.'
        matPath = fullfile(SAMPLES_DIR, items(k).name, 'sample.mat');
        if isfile(matPath)
            sampleDirs{end+1, 1} = items(k).name; %#ok<AGROW>
        end
    end
end

n = numel(sampleDirs);
if n == 0
    error('step9:no_graphs', 'No sample.mat files found under %s', SAMPLES_DIR);
end
fprintf('Found %d graphs. Aggregating...\n', n);

%% Load and aggregate

sample_ids = sampleDirs;   % cell(n,1) of folder names
X_cell     = cell(n, 1);
ei_cell    = cell(n, 1);
N_vec      = zeros(n, 1);
TR_vec     = nan(n, 1);
ANG_vec    = nan(n, 1);
Dense_cell = cell(n, 1);
has_dense_count = 0;

t0 = tic;
for i = 1:n
    d  = load(fullfile(SAMPLES_DIR, sample_ids{i}, 'sample.mat'));
    gd = d.gnn_data;

    X_cell{i}  = single(gd.x.');          % 4 x N  (x, y, radius, boundary)
    ei_cell{i} = double(gd.edge_index);   % 2 x E
    N_vec(i)   = gd.num_nodes;
    TR_vec(i)  = gd.tr_ratio;
    ANG_vec(i) = gd.ang_deg;
    if isfield(d, 'dense_data') && isfield(d.dense_data, 'raster')
        Dense_cell{i} = uint8(d.dense_data.raster);
        has_dense_count = has_dense_count + 1;
    else
        Dense_cell{i} = rasterize_structural_graph(X_cell{i}, RASTER_SIZE);
    end

    if mod(i, 500) == 0
        fprintf('  %d / %d  (%.1f s)\n', i, n, toc(t0));
    end
end
fprintf('Loaded %d graphs in %.1f s\n', n, toc(t0));

%% Save aggregate

if ~isfolder(TARGETS_DIR), mkdir(TARGETS_DIR); end

created_at     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
schema_version = SCHEMA_VERSION;

outFile = fullfile(TARGETS_DIR, 'graphs_all.mat');
save(outFile, 'sample_ids', 'X_cell', 'ei_cell', 'N_vec', 'TR_vec', 'ANG_vec', 'Dense_cell', ...
    'schema_version', 'created_at', '-v7.3');
fprintf('Saved: %s\n', outFile);

%% Manifest

try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

manifest = struct( ...
    'n_graphs',        n, ...
    'n_dense_native',  has_dense_count, ...
    'raster_size',     RASTER_SIZE, ...
    'schema_version',  schema_version, ...
    'git_sha',         git_sha, ...
    'created_at',      created_at);

manifestFile = fullfile(TARGETS_DIR, 'graphs_all_manifest.json');
fid = fopen(manifestFile, 'w');
fprintf(fid, '%s\n', jsonencode(manifest));
fclose(fid);
fprintf('Manifest: %s\n', manifestFile);
fprintf('Done.\n');

function raster = rasterize_structural_graph(X, gridSize)
xy = double(X(1:2, :).');
if size(X, 1) >= 4
    boundary = logical(X(4, :).');
else
    boundary = false(size(xy, 1), 1);
end

xLo = min(xy(:, 1)); xHi = max(xy(:, 1));
yLo = min(xy(:, 2)); yHi = max(xy(:, 2));

cols = 1 + round((xy(:, 1) - xLo) / max(xHi - xLo, eps) * (gridSize - 1));
rows = 1 + round((xy(:, 2) - yLo) / max(yHi - yLo, eps) * (gridSize - 1));
cols = min(gridSize, max(1, cols));
rows = min(gridSize, max(1, rows));

occupancy = false(gridSize, gridSize);
boundaryMask = false(gridSize, gridSize);
idx = sub2ind([gridSize, gridSize], rows, cols);
occupancy(idx) = true;
boundaryMask(idx(boundary)) = true;
raster = uint8(cat(3, occupancy, boundaryMask));
end
