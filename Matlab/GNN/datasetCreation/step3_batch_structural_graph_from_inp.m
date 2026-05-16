% step3_batch_structural_graph_from_inp
%
% Build structural graphs directly from shell INP files.
% This is pre-deformation graph generation (no Abaqus CSV required).
% Output is a single MATLAB sample file per run:
%   Matlab/GNN/data/dataset_hybrid/samples/<run_name>/sample.mat
%
% Input root:
%   Matlab/GNN/data/raw/samples
% Output root:
%   Matlab/GNN/data/dataset_hybrid/samples

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
gnnRoot = fileparts(scriptDir);
repoRoot = fileparts(fileparts(gnnRoot));
helpersDir = fullfile(gnnRoot, 'helpers');

addpath(scriptDir);
addpath(helpersDir);

samplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');

% Pick which detail level to build.
DETAIL_LEVEL = 0.3;

% Pick dataset schema version. 'V1' = legacy 4-dim node features and binary
% edges; 'V2' adds a 5-dim edge_attr and replaces the binary boundary flag
% with a 4-class one-hot {interior, clamped_left, loaded_right, free_top_bot};
% 'V3' adds raster-based boundary node augmentation (radius=0 surface nodes,
% one edge to the nearest centerline node) and replaces boundary-distance
% features with a single dist_to_loaded scalar (node F becomes 8).
% Each version writes to its own folder so older datasets stay intact.
DATASET_VERSION = 'V3';

dsBase = sprintf('dataset_hybrid_d%03d', round(DETAIL_LEVEL * 100));
switch upper(DATASET_VERSION)
    case 'V2', dsBase = [dsBase '_v2'];
    case 'V3', dsBase = [dsBase '_v3'];
end
datasetRoot = fullfile(gnnRoot, 'data', dsBase, 'samples');

ensure_dir(datasetRoot);

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

% GNN-specific minimal graph functions (override originals for this pipeline).
addpath(fullfile(scriptDir, 'gnn_graph'));

OVERWRITE = false;  % skip already-completed samples
MAX_TASKS = inf;    % set small integer for raster/detail smoke exports

runDirs = find_run_dirs(samplesRoot);
fprintf('Step 3 (INP-only): discovered %d run folders.\n', numel(runDirs));

t0 = tic;
nSkipNoInp = 0;
nSkipDone  = 0;
taskInpPaths = {};
taskRunNames = {};

for i = 1:numel(runDirs)
    runDir = runDirs{i};
    [~, runName] = fileparts(runDir);
    inpPath = fullfile(runDir, 'sheet_shell.inp');

    if ~isfile(inpPath)
        nSkipNoInp = nSkipNoInp + 1;
        continue;
    end

    sampleDir = fullfile(datasetRoot, runName);
    sampleMatPath = fullfile(sampleDir, 'sample.mat');
    if ~OVERWRITE && isfile(sampleMatPath)
        nSkipDone = nSkipDone + 1;
        continue;
    end

    taskInpPaths{end+1,1} = inpPath; %#ok<AGROW>
    taskRunNames{end+1,1} = runName; %#ok<AGROW>
end

if isempty(taskInpPaths)
    dt = toc(t0);
    fprintf('Step 3 (INP-only) complete: ok=0 skip_no_inp=%d skip_done=%d fail=0 (%.2f min)\n', ...
        nSkipNoInp, nSkipDone, dt / 60);
    return;
end

if isfinite(MAX_TASKS) && numel(taskInpPaths) > MAX_TASKS
    taskInpPaths = taskInpPaths(1:MAX_TASKS);
    taskRunNames = taskRunNames(1:MAX_TASKS);
    fprintf('Smoke mode: processing first %d graph tasks only.\n', MAX_TASKS);
end

pool = ensure_full_pool();
fprintf('Using parallel pool with %d workers for graph extraction.\n', pool.NumWorkers);

nTasks = numel(taskInpPaths);
ok = false(1, nTasks);
err = cell(1, nTasks);

parfor i = 1:nTasks
    try
        process_single_run(taskInpPaths{i}, datasetRoot, taskRunNames{i}, DETAIL_LEVEL, DATASET_VERSION);
        ok(i) = true;
    catch ME
        err{i} = sprintf('[%s] %s', taskRunNames{i}, ME.message);
    end
end

nOk = nnz(ok);
nFail = nTasks - nOk;
dt = toc(t0);
fprintf('Step 3 (INP-only) complete: ok=%d skip_no_inp=%d skip_done=%d fail=%d (%.2f min)\n', ...
    nOk, nSkipNoInp, nSkipDone, nFail, dt / 60);
if nFail > 0
    badIdx = find(~ok);
    nShow = min(20, numel(badIdx));
    fprintf(2, 'Showing %d/%d failures:\n', nShow, numel(badIdx));
    for k = 1:nShow
        fprintf(2, '  %s\n', err{badIdx(k)});
    end
end

function process_single_run(inpPath, datasetRoot, runName, detailLevel, datasetVersion)
inpData = read_abaqus_inp(inpPath);
fullGraph = build_full_reference_graph_gnn(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});

structuralGraph = extract_structural_graph(fullGraph, ...
    'DetailLevel', detailLevel, ...
    'MinIslandNodes', 1);

sampleDir = fullfile(datasetRoot, runName);
ensure_dir(sampleDir);

sampleMatPath = fullfile(sampleDir, 'sample.mat');

% Dense raster is needed up-front for V3 (it drives the raster-boundary
% augmentation). For V1/V2 the same raster is also saved alongside the
% sample for the CNN branch.
dense_data = rasterize_full_reference_graph(fullGraph, 128);

% Boundary classification from bounding box. V1 stores a single binary flag,
% V2/V3 store a 4-class one-hot {interior, clamped_left, loaded_right,
% free_top_bot}. Corner nodes are assigned in the listed priority order so
% the encoding stays mutually exclusive.
xy  = structuralGraph.node_coords;   % N×2
rng_xy = max(range(xy), [], 1);
tol    = 1e-4 * max(rng_xy);
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));

onLeft   = xy(:,1) - xLo < tol;
onRight  = xHi - xy(:,1) < tol;
onTopBot = (xy(:,2) - yLo < tol) | (yHi - xy(:,2) < tol);

cls_clamped  = onLeft;
cls_loaded   = onRight  & ~cls_clamped;
cls_free     = onTopBot & ~cls_clamped & ~cls_loaded;
cls_interior = ~(cls_clamped | cls_loaded | cls_free);

% Parse TR ratio and loading angle from folder name
tk = regexp(runName, 'tr(\d+)', 'tokens', 'once');
tr_ratio = NaN; if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end
tk = regexp(runName, 'ang(\d+)', 'tokens', 'once');
ang_deg  = NaN; if ~isempty(tk), ang_deg  = str2double(tk{1}); end

gnn_data = struct();
gnn_data.edge_index = int32(structuralGraph.edge_index);
gnn_data.num_nodes  = structuralGraph.num_nodes;
gnn_data.tr_ratio   = tr_ratio;
gnn_data.ang_deg    = ang_deg;
gnn_data.detail_level = detailLevel;

switch upper(datasetVersion)
case 'V2'
    bcOneHot = zeros(size(xy,1), 4);
    bcOneHot(cls_interior, 1) = 1;
    bcOneHot(cls_clamped,  2) = 1;
    bcOneHot(cls_loaded,   3) = 1;
    bcOneHot(cls_free,     4) = 1;

    gnn_data.x = [structuralGraph.node_coords, structuralGraph.node_radius, bcOneHot];  % N×7

    ei = double(structuralGraph.edge_index);   % 2×E
    edge_attr = build_v2_edge_attr(ei, xy, structuralGraph.node_radius(:));
    gnn_data.edge_attr            = single(edge_attr);
    gnn_data.schema_version       = 4;
    gnn_data.schema_version_label = 'V2';

case 'V3'
    % --- Per-wall boundary augmentation. For each wall, group occupied
    %     perimeter pixels by their nearest centerline node:
    %       k >= 2 (perpendicular contact / fan-out): add ONE boundary
    %               node at the group centroid with a single spur to c.
    %       k == 1 (parallel contact / no fan-out): relabel c directly
    %               with the wall's class instead of adding a node.
    %     Centerline nodes that aren't the nearest of any wall pixel keep
    %     class 'interior'. Relabel priority: clamp(2) > load(3) > free(4).
    [newXY, newCls, newParent, relIdx, relCls] = ...
        compute_v3_boundary_augmentation(xy, dense_data);

    classIdxCenter = ones(size(xy, 1), 1);   % interior by default
    for ri = 1:numel(relIdx)
        cIdx  = relIdx(ri);
        newC  = relCls(ri);
        if classIdxCenter(cIdx) == 1 || newC < classIdxCenter(cIdx)
            classIdxCenter(cIdx) = newC;
        end
    end

    Nb = size(newXY, 1);
    Nc = size(xy, 1);

    radius_c = structuralGraph.node_radius(:);
    xy_aug     = [xy; newXY];
    radius_aug = [radius_c; zeros(Nb, 1)];

    % --- New edges: each new boundary node b connects to its parent c.
    %     Stored as (c; b) so (rn - rm)/L = (rc - rb)/L, matching the spec's
    %     unit_dir direction for the new boundary edges.
    ei_old = double(structuralGraph.edge_index);   % 2×E_old
    newEdges = zeros(2, Nb);
    for b = 1:Nb
        newEdges(:, b) = [newParent(b); Nc + b];
    end
    ei_aug = [ei_old, newEdges];

    edge_attr_aug = build_v2_edge_attr(ei_aug, xy_aug, radius_aug);

    classIdxAll          = [classIdxCenter; newCls];
    bcOneHot             = zeros(Nc + Nb, 4);
    linIdx               = sub2ind([Nc + Nb, 4], (1:Nc+Nb)', classIdxAll);
    bcOneHot(linIdx)     = 1;

    xExtent = max(xHi - xLo, eps);
    dist_to_loaded = (xHi - xy_aug(:, 1)) / xExtent;   % 0 at right, 1 at left

    gnn_data.x                    = [xy_aug, radius_aug, bcOneHot, dist_to_loaded];  % N×8
    gnn_data.edge_index           = int32(ei_aug);
    gnn_data.num_nodes            = Nc + Nb;
    gnn_data.edge_attr            = single(edge_attr_aug);
    gnn_data.schema_version       = 5;
    gnn_data.schema_version_label = 'V3';

otherwise   % 'V1'
    boundaryFlag = double(~cls_interior);
    gnn_data.x = [structuralGraph.node_coords, structuralGraph.node_radius, boundaryFlag];  % N×4
    gnn_data.schema_version       = 3;
    gnn_data.schema_version_label = 'V1';
end

save(sampleMatPath, 'gnn_data', 'dense_data', '-v7');
end

function edge_attr = build_v2_edge_attr(ei, xy, radius)
% Per-edge 5-dim feature: [L, unit_dir_x, unit_dir_y, mean_r, |dr|], using
% the (rn - rm) convention where ei(1,:) = n, ei(2,:) = m.
E = size(ei, 2);
if E == 0
    edge_attr = zeros(0, 5);
    return;
end
rn   = xy(ei(1,:), :);
rm   = xy(ei(2,:), :);
diff = rn - rm;
L    = sqrt(sum(diff.^2, 2));
Lsf  = max(L, 1e-12);
unitDir = diff ./ Lsf;
rn_r = radius(ei(1,:));
rm_r = radius(ei(2,:));
mean_r = (rn_r + rm_r) / 2;
abs_dr = abs(rn_r - rm_r);
edge_attr = [L, unitDir, mean_r, abs_dr];   % E×5
end

function [newXY, newCls, newParent, relabelIdx, relabelCls] = ...
        compute_v3_boundary_augmentation(xy, dense_data)
% Per-wall raster-boundary augmentation.
%   For each wall (left/clamp, right/load, bot/free, top/free):
%     - gather occupied perimeter pixels on that wall
%     - find each pixel's nearest centerline node
%     - group pixels by parent node
%     - k>=2 groups (perpendicular contact): emit ONE new boundary node
%       at the group centroid with a single edge to the parent
%     - k==1 groups (parallel contact): no new node, queue the parent
%       for relabeling with the wall's class
% Returns:
%   newXY     Nnew x 2  new boundary node xy
%   newCls    Nnew x 1  class code per new node {2=clamp, 3=load, 4=free}
%   newParent Nnew x 1  centerline node each new node attaches to
%   relabelIdx, relabelCls  pairs of (c, new class) to apply with priority
%                           (clamp 2 < load 3 < free 4).
occ = logical(dense_data.raster(:, :, 1));
[H, W] = size(occ);
xLo = dense_data.xy_bounds(1); xHi = dense_data.xy_bounds(2);
yLo = dense_data.xy_bounds(3); yHi = dense_data.xy_bounds(4);
dxPix = (xHi - xLo) / max(W - 1, 1);
dyPix = (yHi - yLo) / max(H - 1, 1);
relabelMaxDist = 1.5 * max(dxPix, dyPix);   % only relabel c if it sits on the wall

walls = struct('xy', {}, 'cls', {});
rL = find(occ(:, 1));
if ~isempty(rL)
    walls(end+1).xy = [repmat(xLo, numel(rL), 1), yLo + (rL - 1) * dyPix];
    walls(end).cls  = 2;
end
rR = find(occ(:, W));
if ~isempty(rR)
    walls(end+1).xy = [repmat(xHi, numel(rR), 1), yLo + (rR - 1) * dyPix];
    walls(end).cls  = 3;
end
cB = find(occ(1, :));
if ~isempty(cB)
    walls(end+1).xy = [xLo + (cB(:) - 1) * dxPix, repmat(yLo, numel(cB), 1)];
    walls(end).cls  = 4;
end
cT = find(occ(H, :));
if ~isempty(cT)
    walls(end+1).xy = [xLo + (cT(:) - 1) * dxPix, repmat(yHi, numel(cT), 1)];
    walls(end).cls  = 4;
end

newXY      = zeros(0, 2);
newCls     = zeros(0, 1);
newParent  = zeros(0, 1);
relabelIdx = zeros(0, 1);
relabelCls = zeros(0, 1);

for w = 1:numel(walls)
    pxy = walls(w).xy;
    cls = walls(w).cls;
    k   = size(pxy, 1);
    parent = zeros(k, 1);
    for p = 1:k
        d2 = sum((xy - pxy(p, :)).^2, 2);
        [~, parent(p)] = min(d2);
    end
    [uParent, ~, grp] = unique(parent);
    for g = 1:numel(uParent)
        sel = find(grp == g);
        cIdx = uParent(g);
        if numel(sel) >= 2
            newXY(end+1, :)     = mean(pxy(sel, :), 1); %#ok<AGROW>
            newCls(end+1, 1)    = cls;                  %#ok<AGROW>
            newParent(end+1, 1) = cIdx;                 %#ok<AGROW>
        else
            % Singleton: only treat as parallel contact (relabel c) if c is
            % actually on the wall. Otherwise emit a new boundary node at
            % the pixel so we don't drag an interior class into the graph.
            pPos = pxy(sel, :);
            if norm(xy(cIdx, :) - pPos) <= relabelMaxDist
                relabelIdx(end+1, 1) = cIdx;            %#ok<AGROW>
                relabelCls(end+1, 1) = cls;             %#ok<AGROW>
            else
                newXY(end+1, :)     = pPos;             %#ok<AGROW>
                newCls(end+1, 1)    = cls;              %#ok<AGROW>
                newParent(end+1, 1) = cIdx;             %#ok<AGROW>
            end
        end
    end
end
end

function dense_data = rasterize_full_reference_graph(fullGraph, gridSize)
xy = double(fullGraph.node_coords(:, 1:2));
boundary = logical(fullGraph.boundary_mask(:));

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

if isfield(fullGraph, 'edges_local') && ~isempty(fullGraph.edges_local)
    edges = double(fullGraph.edges_local);
    for e = 1:size(edges, 1)
        s = edges(e, 1);
        t = edges(e, 2);
        nSteps = max(abs(rows(s) - rows(t)), abs(cols(s) - cols(t))) + 1;
        rr = round(linspace(rows(s), rows(t), nSteps));
        cc = round(linspace(cols(s), cols(t), nSteps));
        edgeIdx = sub2ind([gridSize, gridSize], rr, cc);
        occupancy(edgeIdx) = true;
    end
end

dense_data = struct();
dense_data.raster = uint8(cat(3, occupancy, boundaryMask));
dense_data.channels = {'occupancy', 'boundary'};
dense_data.grid_size = gridSize;
dense_data.xy_bounds = [xLo, xHi, yLo, yHi];
dense_data.schema_version = 1;
end

function ensure_dir(pathStr)
if ~isfolder(pathStr)
    mkdir(pathStr);
end
end


function out = find_run_dirs(rootDir)
out = {};
if ~isfolder(rootDir)
    return;
end
items = dir(rootDir);
for k = 1:numel(items)
    if ~items(k).isdir || items(k).name(1) == '.'
        continue;
    end
    sub = fullfile(rootDir, items(k).name);
    if isfile(fullfile(sub, 'sheet.mat')) && isfile(fullfile(sub, 'mesh_manifest.json'))
        out{end+1,1} = sub; %#ok<AGROW>
    else
        out = [out; find_run_dirs(sub)]; %#ok<AGROW>
    end
end
end
