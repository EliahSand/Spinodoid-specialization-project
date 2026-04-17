% step1_generate_psscone_shell_dataset
%
%
%
% Output folders are written under:
%   Matlab/GNN/data/raw/samples/trXX/angYYY/<run-folder>

scriptPath = mfilename('fullpath');
scriptDir = fileparts(scriptPath);
gnnRoot = fileparts(scriptDir);
repoRoot = fileparts(fileparts(gnnRoot));
helpersDir = fullfile(gnnRoot, 'helpers');

addpath(scriptDir);
addpath(helpersDir);
addpath(fullfile(repoRoot, 'Matlab'));

% ---------------- Configuration ----------------
numSamples = 6000;
anglesDeg = 0:90;
%anglesDeg = [0 30 60 90];
ratios = [0.5, 1.0, 2.0];
%ratios = [1.0];

baseThickness = 2e-3;                 % m
outRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');

% PSSCone defaults used for dataset generation
baseParams = struct();
baseParams.N = 128;
baseParams.L = 40e-3;
baseParams.lambda_vox = 25;
baseParams.bandwidth = 2;
baseParams.nModes = 4000;
baseParams.solid_frac = 0.50;
baseParams.coneDeg = [30 0 0];
baseParams.slice_count = 8;
baseParams.tilesXY = [1 1];
baseParams.remove_top_spin_frac = 0.0;
% ----------------------------------------------

if ~isfolder(outRoot)
    mkdir(outRoot);
end

nCombos = numel(anglesDeg) * numel(ratios);
if numSamples < nCombos
    error('numSamples (%d) must be >= number of angle-ratio combos (%d).', numSamples, nCombos);
end

basePerCombo = floor(numSamples / nCombos);
extraCombos = mod(numSamples, nCombos);

fprintf('Step 1: target %d samples | combos=%d | base=%d | extra=%d\n', ...
    numSamples, nCombos, basePerCombo, extraCombos);

comboIdx = 0;
t0 = tic;

% Precompute all jobs first, then run one parfor over all jobs.
taskRatio = zeros(1, numSamples);
taskAngle = zeros(1, numSamples);
taskSeed = zeros(1, numSamples);
taskOutRoot = cell(1, numSamples);
taskTrLabel = cell(1, numSamples);
taskAngLabel = cell(1, numSamples);
taskIdx = 0;

for rIdx = 1:numel(ratios)
    ratio = ratios(rIdx);
    trLabel = sprintf('tr%02d', round(100 * ratio));

    for aIdx = 1:numel(anglesDeg)
        ang = anglesDeg(aIdx);
        comboIdx = comboIdx + 1;

        nSeeds = basePerCombo + double(comboIdx <= extraCombos);
        angLabel = sprintf('ang%03d', round(ang));
        comboRoot = fullfile(outRoot, trLabel, angLabel);
        if ~isfolder(comboRoot)
            mkdir(comboRoot);
        end

        fprintf('[%s | %s] generating %d seeds...\n', trLabel, angLabel, nSeeds);
        for seed = 1:nSeeds
            taskIdx = taskIdx + 1;
            taskRatio(taskIdx) = ratio;
            taskAngle(taskIdx) = ang;
            taskSeed(taskIdx) = seed;             % unique per (angle, ratio)
            taskOutRoot{taskIdx} = comboRoot;
            taskTrLabel{taskIdx} = trLabel;
            taskAngLabel{taskIdx} = angLabel;
        end
    end
end

if taskIdx ~= numSamples
    error('Internal error: expected %d tasks but built %d.', numSamples, taskIdx);
end

pool = ensure_full_pool();
fprintf('Using parallel pool with %d workers.\n', pool.NumWorkers);

nTasks = taskIdx;
ok = false(1, nTasks);
err = cell(1, nTasks);

parfor i = 1:nTasks
    p = baseParams;
    p.t_base = baseThickness;
    p.t_spin = taskRatio(i) * baseThickness;
    p.lamellarAngleDeg = taskAngle(i);
    p.rngSeed = taskSeed(i);
    p.resultsRoot = taskOutRoot{i};

    try
        PSSCone(p);
        ok(i) = true;
    catch ME
        err{i} = sprintf('[%s | %s | seed=%d] %s', ...
            taskTrLabel{i}, taskAngLabel{i}, taskSeed(i), ME.message);
    end
end

generated = nnz(ok);
nFail = nTasks - generated;
dt = toc(t0);
fprintf('Step 1 complete: generated=%d failed=%d in %.2f min.\n', generated, nFail, dt / 60);

if nFail > 0
    badIdx = find(~ok);
    nShow = min(20, numel(badIdx));
    fprintf(2, 'Showing %d/%d failures:\n', nShow, numel(badIdx));
    for k = 1:nShow
        fprintf(2, '  %s\n', err{badIdx(k)});
    end
end
