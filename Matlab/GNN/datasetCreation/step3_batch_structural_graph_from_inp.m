% step3_batch_structural_graph_from_inp
%
% Build structural graphs directly from shell INP files.
% This is pre-deformation graph generation (no Abaqus CSV required).
% Output is a single MATLAB sample file per run:
%   Matlab/GNN/data/dataset/samples/<run_name>/sample.mat
%
% Input root:
%   Matlab/GNN/data/raw/samples
% Output root:
%   Matlab/GNN/data/dataset/samples

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
gnnRoot = fileparts(scriptDir);
repoRoot = fileparts(fileparts(gnnRoot));
helpersDir = fullfile(gnnRoot, 'helpers');

addpath(scriptDir);
addpath(helpersDir);

samplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');
datasetRoot = fullfile(gnnRoot, 'data', 'dataset', 'samples');
ensure_dir(datasetRoot);

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

% GNN-specific minimal graph functions (override originals for this pipeline).
addpath(fullfile(scriptDir, 'gnn_graph'));

OVERWRITE = true;   % set false to skip already-completed samples

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
    'DetailLevel', 0.30, ...
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

save(sampleMatPath, 'gnn_data', '-v7');
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
