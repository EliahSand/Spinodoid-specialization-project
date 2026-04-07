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

samplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');
datasetRoot = fullfile(gnnRoot, 'data', 'dataset', 'samples');
ensure_dir(datasetRoot);

gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

runDirs = find_run_dirs(samplesRoot);
fprintf('Step 3 (INP-only): discovered %d run folders.\n', numel(runDirs));

t0 = tic;
nSkip = 0;
taskInpPaths = {};
taskRunNames = {};

for i = 1:numel(runDirs)
    runDir = runDirs{i};
    [~, runName] = fileparts(runDir);
    inpPath = fullfile(runDir, 'sheet_shell.inp');

    if ~isfile(inpPath)
        nSkip = nSkip + 1;
        continue;
    end

    sampleDir = fullfile(datasetRoot, runName);
    sampleMatPath = fullfile(sampleDir, 'sample.mat');
    if isfile(sampleMatPath)
        nSkip = nSkip + 1;
        continue;
    end

    taskInpPaths{end+1,1} = inpPath; %#ok<AGROW>
    taskRunNames{end+1,1} = runName; %#ok<AGROW>
end

if isempty(taskInpPaths)
    dt = toc(t0);
    fprintf('Step 3 (INP-only) complete: ok=0 skip=%d fail=0 (%.2f min)\n', nSkip, dt / 60);
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
fprintf('Step 3 (INP-only) complete: ok=%d skip=%d fail=%d (%.2f min)\n', ...
    nOk, nSkip, nFail, dt / 60);
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
fullGraph = build_full_reference_graph(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});

[structuralGraph, skeletonGraph, ~, ~] = extract_structural_graph(fullGraph, ...
    'DetailLevel', 0.30, ...
    'MinIslandNodes', 1);

[ratioVal, ratioLog2, ratioLabel] = parse_spin_base_ratio(runName);

sampleDir = fullfile(datasetRoot, runName);
ensure_dir(sampleDir);

sampleMatPath = fullfile(sampleDir, 'sample.mat');

sample = struct();
sample.id = runName;
sample.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
sample.source = struct('inp_path', inpPath);
sample.graph = structuralGraph;
sample.graph_info = struct( ...
    'mode', 'structural_graph_predeformation', ...
    'detail_level', 0.30, ...
    'min_island_nodes', 1, ...
    'spin_base_ratio', ratioVal, ...
    'spin_base_ratio_log2', ratioLog2, ...
    'spin_base_ratio_label', ratioLabel, ...
    'full_nodes', fullGraph.num_nodes, ...
    'full_edges', size(fullGraph.edges_local, 1), ...
    'skeleton_nodes', skeletonGraph.num_nodes, ...
    'skeleton_edges', size(skeletonGraph.edges_local, 1), ...
    'structural_nodes', structuralGraph.num_nodes, ...
    'structural_edges', size(structuralGraph.edges_local, 1));

save(sampleMatPath, 'sample', '-v7');
end

function ensure_dir(pathStr)
if ~isfolder(pathStr)
    mkdir(pathStr);
end
end

function [ratioVal, ratioLog2, ratioLabel] = parse_spin_base_ratio(runName)
ratioVal = NaN;
ratioLog2 = NaN;
ratioLabel = '';

tok = regexp(lower(runName), 'tr(\d+)', 'tokens', 'once');
if isempty(tok)
    return;
end

ratioDigits = str2double(tok{1});
if ~isfinite(ratioDigits)
    return;
end

ratioVal = ratioDigits / 100;
if ratioVal > 0
    ratioLog2 = log2(ratioVal);
end
ratioLabel = sprintf('tr%02d', round(ratioDigits));
end

function out = find_run_dirs(rootDir)
out = {};
if ~isfolder(rootDir)
    return;
end

allPath = genpath(rootDir);
dirs = strsplit(allPath, pathsep);
dirs = dirs(~cellfun('isempty', dirs));
for k = 1:numel(dirs)
    d = dirs{k};
    if isfile(fullfile(d, 'sheet.mat')) && isfile(fullfile(d, 'mesh_manifest.json'))
        out{end+1,1} = d; %#ok<AGROW>
    end
end

if ~isempty(out)
    out = unique(out, 'stable');
end
end

function pool = ensure_full_pool()
c = parcluster('local');
targetWorkers = c.NumWorkers;
pool = gcp('nocreate');

if ~isempty(pool) && pool.NumWorkers == targetWorkers
    return;
end

if ~isempty(pool)
    delete(pool);
end

try
    pool = parpool(c, targetWorkers);
catch ME
    warning('Could not start pool with %d workers (%s). Falling back to default pool size.', ...
        targetWorkers, ME.message);
    pool = parpool(c);
end
end
