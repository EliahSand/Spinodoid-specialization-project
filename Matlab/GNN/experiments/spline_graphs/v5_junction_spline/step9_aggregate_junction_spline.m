% step9_aggregate_junction_spline.m
%
% Aggregates schema v5 junction-spline graphs into a single graphs_all.mat.
% Writes under dataset_junction_spline; leaves v3 and v4 datasets untouched.

if ~exist('REPO_ROOT', 'var')
    REPO_ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
if ~exist('SAMPLES_DIR', 'var')
    SAMPLES_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_junction_spline', 'samples');
end
if ~exist('TARGETS_DIR', 'var')
    TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset_junction_spline', 'targets');
end

SCHEMA_VERSION = 5;
RASTER_SIZE    = 128;

%% Discover samples

items = dir(SAMPLES_DIR);
sampleDirs = {};
for k = 1:numel(items)
    if items(k).isdir && items(k).name(1) ~= '.'
        mp = fullfile(SAMPLES_DIR, items(k).name, 'sample.mat');
        if isfile(mp), sampleDirs{end+1,1} = items(k).name; end %#ok<AGROW>
    end
end
n = numel(sampleDirs);
if n == 0
    error('step9_jl:no_graphs', 'No sample.mat files found under %s', SAMPLES_DIR);
end
fprintf('Found %d v5 junction-spline graphs. Aggregating...\n', n);

%% Load + aggregate

sample_ids       = sampleDirs;
X_cell           = cell(n, 1);
EdgeAttr_cell    = cell(n, 1);
ei_cell          = cell(n, 1);
N_vec            = zeros(n, 1);
E_vec            = zeros(n, 1);
TR_vec           = nan(n, 1);
ANG_vec          = nan(n, 1);
Dense_cell       = cell(n, 1);
feature_names      = {};
edge_feature_names = {};
representation   = 'junction_spline_graph';
has_dense_count  = 0;

t0 = tic;
for i = 1:n
    d  = load(fullfile(SAMPLES_DIR, sample_ids{i}, 'sample.mat'));
    gd = d.gnn_data;

    if ~isfield(gd, 'schema_version') || gd.schema_version ~= SCHEMA_VERSION
        error('step9_jl:bad_schema', 'Expected schema_version=%d in %s.', SCHEMA_VERSION, sample_ids{i});
    end

    X_cell{i}        = single(gd.x.');          % F x N
    ei_cell{i}       = double(gd.edge_index);   % 2 x E
    EdgeAttr_cell{i} = single(gd.edge_features.');  % Fe x E
    N_vec(i)         = gd.num_nodes;
    E_vec(i)         = size(gd.edge_features, 1);
    TR_vec(i)        = gd.tr_ratio;
    ANG_vec(i)       = gd.ang_deg;

    if isempty(feature_names) && isfield(gd, 'feature_names')
        feature_names = gd.feature_names;
    end
    if isempty(edge_feature_names) && isfield(gd, 'edge_feature_names')
        edge_feature_names = gd.edge_feature_names;
    end

    if isfield(d, 'dense_data') && isfield(d.dense_data, 'raster')
        Dense_cell{i} = uint8(d.dense_data.raster);
        has_dense_count = has_dense_count + 1;
    else
        Dense_cell{i} = rasterize_xy(X_cell{i}, RASTER_SIZE);
    end

    if mod(i, 500) == 0
        fprintf('  %d / %d  (%.1f s)\n', i, n, toc(t0));
    end
end
fprintf('Loaded %d graphs in %.1f s\n', n, toc(t0));

%% Save

if ~isfolder(TARGETS_DIR), mkdir(TARGETS_DIR); end
created_at     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
schema_version = SCHEMA_VERSION;
outFile = fullfile(TARGETS_DIR, 'graphs_all.mat');

save(outFile, 'sample_ids', 'X_cell', 'EdgeAttr_cell', 'ei_cell', ...
    'N_vec', 'E_vec', 'TR_vec', 'ANG_vec', 'Dense_cell', ...
    'feature_names', 'edge_feature_names', 'representation', ...
    'schema_version', 'created_at', '-v7.3');
fprintf('Saved: %s\n', outFile);

%% Manifest

try, [~, sha] = system('git rev-parse HEAD'); sha = strtrim(sha);
catch, sha = 'unknown'; end

manifest = struct('n_graphs', n, 'n_dense_native', has_dense_count, ...
    'raster_size', RASTER_SIZE, 'schema_version', schema_version, ...
    'representation', representation, 'feature_names', {feature_names}, ...
    'edge_feature_names', {edge_feature_names}, 'git_sha', sha, 'created_at', created_at);
mf = fullfile(TARGETS_DIR, 'graphs_all_manifest.json');
fid = fopen(mf, 'w'); fprintf(fid, '%s\n', jsonencode(manifest)); fclose(fid);
fprintf('Manifest: %s\nDone.\n', mf);

function raster = rasterize_xy(X, gridSize)
xy = double(X(1:2,:).');
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));
cols = 1 + round((xy(:,1)-xLo)/max(xHi-xLo,eps)*(gridSize-1));
rows = 1 + round((xy(:,2)-yLo)/max(yHi-yLo,eps)*(gridSize-1));
cols = min(gridSize,max(1,cols)); rows = min(gridSize,max(1,rows));
occ  = false(gridSize,gridSize);
occ(sub2ind([gridSize,gridSize],rows,cols)) = true;
raster = uint8(cat(3, occ, false(gridSize,gridSize)));
end
