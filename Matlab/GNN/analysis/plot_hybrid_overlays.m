%% Plot true-vs-predicted overlays for a completed GNN run.
%
% Works for any modelMode saved by Main_hybrid_vision.m
% ('hybrid', 'dense_only', 'graph_only'). Uses saved predictions from
% final_metrics.mat. It does not rerun inference or retrain the model.

clear; close all; clc;

%% Config

SPLIT = 'train';          % 'train', 'val', or 'test'
N_SHOW = 12;
SELECTION = 'median';     % 'even', 'first', 'random', 'best', 'median', or 'worst'
RNG_SEED = 2;
MODE = 'dense_only';           % 'hybrid', 'dense_only', 'graph_only', or 'any'
RUN_DIR = '';           % empty = newest completed run matching MODE and current targets

scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fileparts(scriptDir);

targetsDir = fullfile(gnnRoot, 'data', 'dataset', 'targets');

%% Load run and target data

pt = load(fullfile(targetsDir, 'pca_targets.mat'), ...
    'Z_train', 'Z_val', 'Z_test', 'train_idx', 'val_idx', 'test_idx');
pm = load(fullfile(targetsDir, 'pca_model.mat'), 's_grid');
ut = load(fullfile(targetsDir, 'u3_targets.mat'), 'U3_mat', 'sample_ids');

if isempty(RUN_DIR)
    RUN_DIR = newest_compatible_run(fullfile(gnnRoot, 'models', 'best'), pt, MODE);
end

metricsFile = fullfile(RUN_DIR, 'final_metrics.mat');
if ~isfile(metricsFile)
    error('plot_hybrid_overlays:missingMetrics', ...
        'final_metrics.mat not found: %s', metricsFile);
end

mf = load(metricsFile);

if isfield(mf, 'cfg') && isfield(mf.cfg, 'modelMode')
    runMode = mf.cfg.modelMode;
else
    runMode = infer_mode_from_dirname(RUN_DIR);
end

hasSavedIdx = all(isfield(mf, {'train_idx', 'val_idx', 'test_idx'}));
isSubset = isfield(mf, 'cfg') && isfield(mf.cfg, 'subsetFraction') && mf.cfg.subsetFraction < 1;
if isSubset && ~hasSavedIdx
    error('plot_hybrid_overlays:subsetRun', ...
        ['Run has subsetFraction=%.3g but final_metrics.mat does not contain ', ...
         'train_idx/val_idx/test_idx. Re-train with the updated Main_hybrid_vision.m.'], ...
        mf.cfg.subsetFraction);
end

[metrics, Z_true, split_idx, splitLabel] = select_split(SPLIT, mf, pt);
U3_true = ut.U3_mat(split_idx, :);
U3_pred = metrics.U3_pred;
Z_pred = metrics.Zhat;

check_run_matches_targets(U3_true, U3_pred, Z_true, Z_pred, RUN_DIR, splitLabel);

errU3 = U3_pred - U3_true;
u3Rmse = sqrt(mean(errU3.^2, 2));
zRmse = sqrt(mean((Z_pred - Z_true).^2, 2));

showIdx = choose_samples(SELECTION, N_SHOW, u3Rmse, RNG_SEED);
sampleIds = ut.sample_ids(split_idx);

outDir = fullfile(RUN_DIR, 'analysis', splitLabel);
if ~isfolder(outDir), mkdir(outDir); end

%% Export overlays

modeLabel = mode_display_label(runMode);
splitTag = sprintf('%s [%s]', splitLabel, modeLabel);

u3File = fullfile(outDir, 'u3_true_vs_pred_overlay.png');
plot_u3_overlay(u3File, pm.s_grid, U3_true, U3_pred, u3Rmse, sampleIds, showIdx, splitTag);

pcaOverlayFile = fullfile(outDir, 'pca_true_vs_pred_overlay.png');
plot_pca_overlay(pcaOverlayFile, Z_true, Z_pred, zRmse, sampleIds, showIdx, splitTag);

pcaScatterFile = fullfile(outDir, 'pca_component_scatter.png');
plot_pca_scatter(pcaScatterFile, Z_true, Z_pred, splitTag);

fprintf('%s overlay plots written to:\n  %s\n', runMode, outDir);
fprintf('  %s\n', u3File);
fprintf('  %s\n', pcaOverlayFile);
fprintf('  %s\n', pcaScatterFile);

%% Local functions

function [metrics, Z_true, split_idx, splitLabel] = select_split(splitName, mf, pt)
splitLabel = lower(string(splitName));
switch splitLabel
    case "train"
        metrics = mf.trainMetrics;
        [Z_true, split_idx] = resolve_truth(mf, pt, 'train');
    case "val"
        metrics = mf.valMetrics;
        [Z_true, split_idx] = resolve_truth(mf, pt, 'val');
    case "test"
        metrics = mf.testMetrics;
        [Z_true, split_idx] = resolve_truth(mf, pt, 'test');
    otherwise
        error('plot_hybrid_overlays:badSplit', ...
            'SPLIT must be train, val, or test. Got: %s', splitName);
end
splitLabel = char(splitLabel);
end

function [Z_true, split_idx] = resolve_truth(mf, pt, splitName)
ptIdxField = sprintf('%s_idx', splitName);
ptZField = sprintf('Z_%s', splitName);
ptFull = pt.(ptIdxField);
ZFull = pt.(ptZField);

if isfield(mf, ptIdxField)
    split_idx = mf.(ptIdxField);
    [tf, pos] = ismember(split_idx, ptFull);
    if ~all(tf)
        error('plot_hybrid_overlays:idxMismatch', ...
            ['Saved %s indices contain entries not in current pca_targets.mat. ', ...
             'The dataset split has changed since this run was trained.'], ptIdxField);
    end
    Z_true = ZFull(pos, :);
else
    split_idx = ptFull;
    Z_true = ZFull;
end
end

function runDir = newest_compatible_run(bestRoot, pt, mode)
modeLower = lower(string(mode));
allowedModes = ["hybrid", "dense_only", "graph_only", "any"];
if ~any(modeLower == allowedModes)
    error('plot_hybrid_overlays:badMode', ...
        'MODE must be one of hybrid, dense_only, graph_only, any. Got: %s', mode);
end

if modeLower == "any"
    pattern = '*';
else
    pattern = sprintf('*_%s', char(modeLower));
end

dirs = dir(fullfile(bestRoot, pattern));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
if isempty(dirs)
    error('plot_hybrid_overlays:noRuns', ...
        'No run directories matching "%s" found under %s', pattern, bestRoot);
end

targetCounts = [size(pt.Z_train, 1), size(pt.Z_val, 1), size(pt.Z_test, 1)];
[~, order] = sort([dirs.datenum], 'descend');

for ii = order
    candidate = fullfile(bestRoot, dirs(ii).name);
    metricsFile = fullfile(candidate, 'final_metrics.mat');
    if ~isfile(metricsFile)
        continue;
    end
    s = load(metricsFile);
    predCounts = [size(s.trainMetrics.U3_pred, 1), ...
                  size(s.valMetrics.U3_pred, 1), ...
                  size(s.testMetrics.U3_pred, 1)];

    if all(isfield(s, {'train_idx', 'val_idx', 'test_idx'}))
        idxOk = all(ismember(s.train_idx, pt.train_idx)) && ...
                all(ismember(s.val_idx,   pt.val_idx))   && ...
                all(ismember(s.test_idx,  pt.test_idx))  && ...
                numel(s.train_idx) == predCounts(1) && ...
                numel(s.val_idx)   == predCounts(2) && ...
                numel(s.test_idx)  == predCounts(3);
        if idxOk
            runDir = candidate;
            fprintf('RUN_DIR empty; using newest compatible run (mode=%s, subset=%d/%d/%d):\n  %s\n', ...
                char(modeLower), predCounts(1), predCounts(2), predCounts(3), runDir);
            return;
        end
    elseif isequal(predCounts, targetCounts)
        runDir = candidate;
        fprintf('RUN_DIR empty; using newest compatible run (mode=%s, full split):\n  %s\n', ...
            char(modeLower), runDir);
        return;
    end
end

error('plot_hybrid_overlays:noCompatibleRun', ...
    ['No completed run matching mode="%s" under %s matches the current target counts ', ...
     '[train=%d val=%d test=%d]. Set RUN_DIR manually for a compatible run.'], ...
    char(modeLower), bestRoot, targetCounts(1), targetCounts(2), targetCounts(3));
end

function modeName = infer_mode_from_dirname(runDir)
[~, name] = fileparts(runDir);
if endsWith(name, '_hybrid')
    modeName = 'hybrid';
elseif endsWith(name, '_dense_only')
    modeName = 'dense_only';
elseif endsWith(name, '_graph_only')
    modeName = 'graph_only';
else
    modeName = 'unknown';
end
end

function label = mode_display_label(modeName)
switch lower(string(modeName))
    case "hybrid",     label = 'hybrid (GNN+CNN)';
    case "dense_only", label = 'CNN only';
    case "graph_only", label = 'GNN only';
    otherwise,         label = char(modeName);
end
end

function check_run_matches_targets(U3_true, U3_pred, Z_true, Z_pred, runDir, splitLabel)
if size(U3_true, 1) ~= size(U3_pred, 1) || size(Z_true, 1) ~= size(Z_pred, 1)
    error('plot_hybrid_overlays:datasetMismatch', ...
        ['Run predictions do not match the current %s targets.\n', ...
         '  RUN_DIR: %s\n', ...
         '  U3 true/pred samples: %d/%d\n', ...
         '  PCA true/pred samples: %d/%d\n', ...
         'This usually means the run was made before the dataset split changed. ', ...
         'Use a run evaluated against the current targets, or restore the matching target files.'], ...
        splitLabel, runDir, size(U3_true, 1), size(U3_pred, 1), ...
        size(Z_true, 1), size(Z_pred, 1));
end
end

function idx = choose_samples(selection, nShow, rmse, rngSeed)
n = numel(rmse);
nShow = min(max(1, nShow), n);
selection = lower(string(selection));

switch selection
    case "even"
        idx = unique(round(linspace(1, n, nShow)), 'stable');
    case "first"
        idx = 1:nShow;
    case "random"
        rng(rngSeed, 'twister');
        idx = sort(randperm(n, nShow));
    case "best"
        [~, order] = sort(rmse, 'ascend');
        idx = order(1:nShow);
    case "median"
        med = median(rmse);
        [~, order] = sort(abs(rmse - med), 'ascend');
        idx = order(1:nShow);
    case "worst"
        [~, order] = sort(rmse, 'descend');
        idx = order(1:nShow);
    otherwise
        error('plot_hybrid_overlays:badSelection', ...
            'SELECTION must be even, first, random, best, median, or worst. Got: %s', selection);
end
idx = idx(:).';
end

function plot_u3_overlay(outFile, s_grid, U3_true, U3_pred, rmse, sampleIds, showIdx, splitLabel)
[nRows, nCols] = grid_shape(numel(showIdx));
fig = figure('Visible', 'off', 'Position', [100 100 1500 900]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for ii = 1:numel(showIdx)
    k = showIdx(ii);
    nexttile;
    plot(s_grid, U3_true(k, :), 'k-', 'LineWidth', 1.2); hold on;
    plot(s_grid, U3_pred(k, :), 'r--', 'LineWidth', 1.2);
    grid on;
    title(sprintf('#%d  %s  RMSE=%.3g', k, compact_id(sampleIds{k}), rmse(k)), ...
        'FontSize', 8, 'Interpreter', 'none');
    xlabel('s');
    ylabel('u3 (m)');
    if ii == 1
        legend('True', 'Pred', 'Location', 'best');
    end
end

sgtitle(sprintf('u3(s) true vs predicted overlays (%s)', splitLabel));
exportgraphics(fig, outFile, 'Resolution', 160);
close(fig);
end

function plot_pca_overlay(outFile, Z_true, Z_pred, rmse, sampleIds, showIdx, splitLabel)
[nRows, nCols] = grid_shape(numel(showIdx));
nComp = size(Z_true, 2);
pc = 1:nComp;

fig = figure('Visible', 'off', 'Position', [100 100 1500 900]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for ii = 1:numel(showIdx)
    k = showIdx(ii);
    nexttile;
    plot(pc, Z_true(k, :), 'ko-', 'LineWidth', 1.0, 'MarkerSize', 3); hold on;
    plot(pc, Z_pred(k, :), 'ro--', 'LineWidth', 1.0, 'MarkerSize', 3);
    grid on;
    xlim([1 nComp]);
    title(sprintf('#%d  %s  RMSE=%.3g', k, compact_id(sampleIds{k}), rmse(k)), ...
        'FontSize', 8, 'Interpreter', 'none');
    xlabel('PC');
    ylabel('coefficient');
    if ii == 1
        legend('True', 'Pred', 'Location', 'best');
    end
end

sgtitle(sprintf('PCA coefficients true vs predicted overlays (%s)', splitLabel));
exportgraphics(fig, outFile, 'Resolution', 160);
close(fig);
end

function plot_pca_scatter(outFile, Z_true, Z_pred, splitLabel)
nComp = size(Z_true, 2);
nRows = 2;
nCols = ceil(nComp / nRows);
fig = figure('Visible', 'off', 'Position', [100 100 1500 700]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:nComp
    nexttile;
    scatter(Z_true(:, c), Z_pred(:, c), 8, 'filled', 'MarkerFaceAlpha', 0.35);
    hold on;
    lo = min([Z_true(:, c); Z_pred(:, c)]);
    hi = max([Z_true(:, c); Z_pred(:, c)]);
    plot([lo hi], [lo hi], 'k--', 'LineWidth', 1.0);
    grid on;
    axis equal tight;
    r2 = component_r2(Z_true(:, c), Z_pred(:, c));
    title(sprintf('PC%d  R2=%.3f', c, r2));
    xlabel('True');
    ylabel('Pred');
end

sgtitle(sprintf('PCA component scatter (%s)', splitLabel));
exportgraphics(fig, outFile, 'Resolution', 160);
close(fig);
end

function r2 = component_r2(yTrue, yPred)
ssRes = sum((yPred - yTrue).^2);
ssTot = sum((yTrue - mean(yTrue)).^2);
r2 = 1 - ssRes / max(ssTot, eps);
end

function [nRows, nCols] = grid_shape(n)
nCols = ceil(sqrt(n));
nRows = ceil(n / nCols);
end

function sid = compact_id(sid)
sid = char(sid);
tk = regexp(sid, 'tr(\d+)_ang(\d+)', 'tokens', 'once');
if ~isempty(tk)
    runTk = regexp(sid, '_run(\d+)', 'tokens', 'once');
    if isempty(runTk)
        sid = sprintf('tr%s ang%s', tk{1}, tk{2});
    else
        sid = sprintf('tr%s ang%s r%s', tk{1}, tk{2}, runTk{1});
    end
    return;
end

maxLen = 24;
if strlength(sid) > maxLen
    sid = ['...' sid(end - maxLen + 4:end)];
end
end
