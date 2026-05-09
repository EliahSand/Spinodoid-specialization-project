% step3_batch_junction_spline
%
% Build junction-spline (v5) graphs from shell INP files.
% Produces one sample.mat per run under dataset_junction_spline/samples/.
% Does not touch dataset_hybrid (v3) or dataset_mat_spline (v4).
%
% Output:
%   Matlab/GNN/data/dataset_junction_spline/samples/<run_name>/sample.mat

scriptPath = mfilename('fullpath');
scriptDir  = fileparts(scriptPath);
gnnRoot    = fileparts(scriptDir);
repoRoot   = fileparts(fileparts(gnnRoot));
helpersDir = fullfile(gnnRoot, 'helpers');

addpath(scriptDir);
addpath(helpersDir);

if ~exist('samplesRoot', 'var')
    samplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');
end
if ~exist('datasetRoot', 'var')
    datasetRoot = fullfile(gnnRoot, 'data', 'dataset_junction_spline', 'samples');
end
ensure_dir(datasetRoot);

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));
addpath(fullfile(scriptDir, 'gnn_graph'));

if ~exist('OVERWRITE', 'var'), OVERWRITE = false; end
if ~exist('MAX_TASKS', 'var'), MAX_TASKS = inf;   end

runDirs = find_run_dirs(samplesRoot);
fprintf('Step 3 junction-spline v5: discovered %d run folders.\n', numel(runDirs));

t0 = tic;
nSkipNoInp = 0; nSkipDone = 0;
taskInpPaths = {}; taskRunNames = {};

for i = 1:numel(runDirs)
    runDir  = runDirs{i};
    [~, runName] = fileparts(runDir);
    inpPath = fullfile(runDir, 'sheet_shell.inp');
    if ~isfile(inpPath), nSkipNoInp = nSkipNoInp + 1; continue; end
    sampleMatPath = fullfile(datasetRoot, runName, 'sample.mat');
    if ~OVERWRITE && isfile(sampleMatPath), nSkipDone = nSkipDone + 1; continue; end
    taskInpPaths{end+1,1} = inpPath;  %#ok<AGROW>
    taskRunNames{end+1,1} = runName;  %#ok<AGROW>
end

if isempty(taskInpPaths)
    dt = toc(t0);
    fprintf('Step 3 v5 complete: ok=0 skip_no_inp=%d skip_done=%d fail=0 (%.2f min)\n', ...
        nSkipNoInp, nSkipDone, dt/60);
    return;
end

if isfinite(MAX_TASKS) && numel(taskInpPaths) > MAX_TASKS
    taskInpPaths = taskInpPaths(1:MAX_TASKS);
    taskRunNames = taskRunNames(1:MAX_TASKS);
    fprintf('Smoke mode: first %d tasks.\n', MAX_TASKS);
end

pool = ensure_full_pool();
fprintf('Using %d parallel workers.\n', pool.NumWorkers);

nTasks = numel(taskInpPaths);
ok  = false(1, nTasks);
err = cell(1, nTasks);

parfor i = 1:nTasks
    try
        process_single_run(taskInpPaths{i}, datasetRoot, taskRunNames{i});
        ok(i) = true;
    catch ME
        err{i} = sprintf('[%s] %s', taskRunNames{i}, ME.message);
    end
end

nOk = nnz(ok); nFail = nTasks - nOk; dt = toc(t0);
fprintf('Step 3 v5 complete: ok=%d skip_no_inp=%d skip_done=%d fail=%d (%.2f min)\n', ...
    nOk, nSkipNoInp, nSkipDone, nFail, dt/60);
if nFail > 0
    badIdx = find(~ok);
    nShow  = min(20, numel(badIdx));
    fprintf(2, 'Showing %d/%d failures:\n', nShow, numel(badIdx));
    for k = 1:nShow, fprintf(2, '  %s\n', err{badIdx(k)}); end
end

% =========================================================================
function process_single_run(inpPath, datasetRoot, runName)
inpData  = read_abaqus_inp(inpPath);
fullGraph = build_full_reference_graph_gnn(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});

[~, ~, dbg] = extract_structural_graph_gnn(fullGraph, ...
    'DetailLevel', 1.00, ...
    'MinIslandNodes', 1);

jlGraph = build_junction_spline_graph( ...
    dbg.skeleton_mask, dbg.radius_map, dbg.boundary_mask, ...
    dbg.x_grid,        dbg.y_grid);

sampleDir = fullfile(datasetRoot, runName);
ensure_dir(sampleDir);

[tr_ratio, ang_deg] = parse_run_name(runName);

gnn_data = struct();
gnn_data.edge_index       = int32(jlGraph.edge_index);
gnn_data.num_nodes        = jlGraph.num_nodes;
gnn_data.x                = jlGraph.node_features;
gnn_data.feature_names    = jlGraph.feature_names;
gnn_data.edge_features    = jlGraph.edge_features;
gnn_data.edge_feature_names = jlGraph.edge_feature_names;
gnn_data.tr_ratio         = tr_ratio;
gnn_data.ang_deg          = ang_deg;
gnn_data.schema_version   = 5;
gnn_data.representation   = 'junction_spline_graph';

dense_data = rasterize_full_reference_graph(fullGraph, 128);

save(fullfile(sampleDir, 'sample.mat'), 'gnn_data', 'dense_data', 'jlGraph', '-v7');
end

function [tr_ratio, ang_deg] = parse_run_name(runName)
tr_ratio = NaN; ang_deg = NaN;
tk = regexp(runName, 'tr(\d+)', 'tokens', 'once');
if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end
tk = regexp(runName, 'ang(\d+)', 'tokens', 'once');
if ~isempty(tk), ang_deg = str2double(tk{1}); end
end

function dense_data = rasterize_full_reference_graph(fullGraph, gridSize)
xy       = double(fullGraph.node_coords(:,1:2));
boundary = logical(fullGraph.boundary_mask(:));
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));
cols = 1 + round((xy(:,1)-xLo)/max(xHi-xLo,eps)*(gridSize-1));
rows = 1 + round((xy(:,2)-yLo)/max(yHi-yLo,eps)*(gridSize-1));
cols = min(gridSize, max(1, cols));
rows = min(gridSize, max(1, rows));
occ  = false(gridSize, gridSize);
bMsk = false(gridSize, gridSize);
idx  = sub2ind([gridSize,gridSize], rows, cols);
occ(idx) = true; bMsk(idx(boundary)) = true;
if isfield(fullGraph,'edges_local') && ~isempty(fullGraph.edges_local)
    edges = double(fullGraph.edges_local);
    for e = 1:size(edges,1)
        s = edges(e,1); t = edges(e,2);
        n = max(abs(rows(s)-rows(t)), abs(cols(s)-cols(t))) + 1;
        rr = round(linspace(rows(s),rows(t),n));
        cc = round(linspace(cols(s),cols(t),n));
        occ(sub2ind([gridSize,gridSize],rr,cc)) = true;
    end
end
dense_data = struct('raster', uint8(cat(3,occ,bMsk)), ...
    'channels', {{'occupancy','boundary'}}, 'grid_size', gridSize, ...
    'xy_bounds', [xLo,xHi,yLo,yHi], 'schema_version', 1);
end

function ensure_dir(p), if ~isfolder(p), mkdir(p); end; end

function out = find_run_dirs(rootDir)
out = {};
if ~isfolder(rootDir), return; end
items = dir(rootDir);
for k = 1:numel(items)
    if ~items(k).isdir || items(k).name(1) == '.', continue; end
    sub = fullfile(rootDir, items(k).name);
    if isfile(fullfile(sub,'sheet.mat')) && isfile(fullfile(sub,'mesh_manifest.json'))
        out{end+1,1} = sub; %#ok<AGROW>
    else
        out = [out; find_run_dirs(sub)]; %#ok<AGROW>
    end
end
end
