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
datasetRoot = fullfile(gnnRoot, 'data', 'dataset_hybrid', 'samples');
ensure_dir(datasetRoot);

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

% GNN-specific minimal graph functions (override originals for this pipeline).
addpath(fullfile(scriptDir, 'gnn_graph'));

OVERWRITE = true;   % set false to skip already-completed samples
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
        process_single_run(taskInpPaths{i}, datasetRoot, taskRunNames{i});
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

function process_single_run(inpPath, datasetRoot, runName)
inpData = read_abaqus_inp(inpPath);
fullGraph = build_full_reference_graph_gnn(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});

structuralGraph = extract_structural_graph_gnn(fullGraph, ...
    'DetailLevel', 1.00, ...
    'MinIslandNodes', 1);

sampleDir = fullfile(datasetRoot, runName);
ensure_dir(sampleDir);

sampleMatPath = fullfile(sampleDir, 'sample.mat');

% Boundary flag: node is on the sheet perimeter if its XY is within
% tolerance of the bounding box. Tolerance = 1e-4 × max(XY range).
xy  = structuralGraph.node_coords;   % N×2
rng_xy = max(range(xy), [], 1);
tol    = 1e-4 * max(rng_xy);
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));
boundaryFlag = double( ...
    xy(:,1) - xLo < tol | xHi - xy(:,1) < tol | ...
    xy(:,2) - yLo < tol | yHi - xy(:,2) < tol);   % N×1

% Parse TR ratio and loading angle from folder name
tk = regexp(runName, 'tr(\d+)', 'tokens', 'once');
tr_ratio = NaN; if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end
tk = regexp(runName, 'ang(\d+)', 'tokens', 'once');
ang_deg  = NaN; if ~isempty(tk), ang_deg  = str2double(tk{1}); end

gnn_data = struct();
gnn_data.edge_index = int32(structuralGraph.edge_index);
gnn_data.num_nodes  = structuralGraph.num_nodes;
gnn_data.x          = [structuralGraph.node_coords, structuralGraph.node_radius, boundaryFlag];  % N×4
gnn_data.tr_ratio   = tr_ratio;
gnn_data.ang_deg    = ang_deg;
gnn_data.detail_level = 1.00;
gnn_data.schema_version = 3;

dense_data = rasterize_full_reference_graph(fullGraph, 128);

save(sampleMatPath, 'gnn_data', 'dense_data', '-v7');
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
