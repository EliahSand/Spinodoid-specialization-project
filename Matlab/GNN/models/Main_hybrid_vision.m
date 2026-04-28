%% Hybrid vision GNN to predict spinodoid sheet deformation PCA targets

clear; close all; clc; clear;
rng(2);

scriptDir = fileparts(mfilename('fullpath'));
gnnRoot = fileparts(scriptDir);
addpath(genpath(fullfile(gnnRoot, 'helpers')));

%% Config

training       = true;
modelMode      = 'hybrid';   % 'hybrid', 'dense_only', or 'graph_only'
subsetFraction = 0.1;        % use <1 only for smoke tests

K              = 4;
hiddenDim      = 96;
cnnChannels    = [32, 64, 128];
fusionDim      = 192;
dropoutRate    = 0.2;
batchSize      = 32;
evalBatchSize  = 32;
maxEpochs      = 180;
lr0            = 7e-4;
decay_rate     = 0.985;
weight_decay   = 1e-5;
valFreq        = 5;
patience       = 12;
minDeltaR2     = 1e-4;
gpuUse         = canUseGPU;

switch modelMode
    case 'hybrid'
        useGraph = true;  useDense = true;
    case 'dense_only'
        useGraph = false; useDense = true;
    case 'graph_only'
        useGraph = true;  useDense = false;
    otherwise
        error('Unknown modelMode: %s', modelMode);
end

runName = sprintf('%s_%s', datestr(now, 'yyyymmdd_HHMMSS'), modelMode);

targetsDir = fullfile(gnnRoot, 'data', 'dataset', 'targets');
samplesDir = fullfile(gnnRoot, 'data', 'dataset_hybrid', 'samples');
hybridTargetsDir = fullfile(gnnRoot, 'data', 'dataset_hybrid', 'targets');
bestDir = fullfile(scriptDir, 'best', runName);
if ~isfolder(bestDir), mkdir(bestDir); end

paramFile = fullfile(bestDir, 'best_params.mat');
progressFile = fullfile(bestDir, 'training_log.mat');

if ~isfile(fullfile(hybridTargetsDir, 'graphs_all.mat'))
    error('Main_hybrid_vision:missingHybridGraphs', ...
        ['Hybrid graphs not found. Run step3_batch_structural_graph_from_inp.m, ', ...
         'then step9_aggregate_graphs.m to create %s.'], ...
        fullfile(hybridTargetsDir, 'graphs_all.mat'));
end

%% Load targets

t = load(fullfile(targetsDir, 'pca_targets.mat'), ...
    'Z_train', 'Z_val', 'Z_test', ...
    'train_ids', 'val_ids', 'test_ids', ...
    'train_idx', 'val_idx', 'test_idx', 'n_components');
p = load(fullfile(targetsDir, 'pca_model.mat'), 'coeff', 'u3_mean', 's_grid', 'explained');
u = load(fullfile(targetsDir, 'u3_targets.mat'), 'U3_mat');

Z_train = t.Z_train;
Z_val = t.Z_val;
Z_test = t.Z_test;
train_ids = t.train_ids;
val_ids = t.val_ids;
test_ids = t.test_ids;
train_idx = t.train_idx;
val_idx = t.val_idx;
test_idx = t.test_idx;
nComp = t.n_components;

coeff = p.coeff;
u3_mean = p.u3_mean;
s_grid = p.s_grid;
explained = p.explained;
U3_mat = u.U3_mat;

%% Optional subset for smoke tests

if subsetFraction < 1.0
    assert(subsetFraction > 0, 'subsetFraction must be > 0');
    nTrSub = max(1, round(numel(train_ids) * subsetFraction));
    nVSub = max(1, round(numel(val_ids) * subsetFraction));
    nTsSub = max(1, round(numel(test_ids) * subsetFraction));

    trSel = sort(randperm(numel(train_ids), nTrSub));
    vSel = sort(randperm(numel(val_ids), nVSub));
    tsSel = sort(randperm(numel(test_ids), nTsSub));

    train_ids = train_ids(trSel); Z_train = Z_train(trSel, :); train_idx = train_idx(trSel);
    val_ids = val_ids(vSel);     Z_val = Z_val(vSel, :);       val_idx = val_idx(vSel);
    test_ids = test_ids(tsSel);  Z_test = Z_test(tsSel, :);    test_idx = test_idx(tsSel);

    fprintf('SUBSET: using %.0f%% -> train=%d val=%d test=%d\n', ...
        100 * subsetFraction, nTrSub, nVSub, nTsSub);
end

nTrain = numel(train_ids);
nVal = numel(val_ids);
nTest = numel(test_ids);

%% Load hybrid graph + raster inputs

all_ids = [train_ids; val_ids; test_ids];
fprintf('Loading %d hybrid samples...\n', numel(all_ids));
[allX, allEI, allN, G_all, Dense_all] = load_hybrid_graph_dataset(all_ids, samplesDir);

nAll = nTrain + nVal + nTest;
trainLocal = (1:nTrain)';
valLocal = (nTrain + 1:nTrain + nVal)';
testLocal = (nTrain + nVal + 1:nAll)';

trainMask = false(nAll, 1);
trainMask(trainLocal) = true;

maxN = max(allN);
fprintf('maxN=%d  nAll=%d  raster=%dx%d\n', maxN, nAll, size(Dense_all, 1), size(Dense_all, 2));

[X_pad, nodeMask, normStats] = pad_and_normalize_hybrid_graphs(allX, allN, trainMask, []);
A_hat = build_norm_adjacency(allEI, allN, maxN);
Dense_input = prepare_dense_raster_inputs(Dense_all);

%% Standardize PCA targets and globals

Z_mean = mean(Z_train, 1);
Z_std = std(Z_train, 0, 1);
Z_std(Z_std < 1e-8) = 1e-8;

normStats.Z_mean = Z_mean;
normStats.Z_std = Z_std;

loss_w = single(sqrt(explained(1:nComp)));
loss_w = loss_w / mean(loss_w);
loss_w = reshape(loss_w, nComp, 1, 1);
normStats.loss_w = loss_w;

Z_all = [Z_train; Z_val; Z_test];
Z_all_std = (Z_all - Z_mean) ./ Z_std;
Z_target = single(reshape(Z_all_std', nComp, 1, []));

tr_mean = mean(G_all(trainLocal, 1));
tr_std = max(std(G_all(trainLocal, 1)), 1e-8);
normStats.tr_mean = tr_mean;
normStats.tr_std = tr_std;

ang_rad = G_all(:, 2) * pi / 180;
tr_std_all = (G_all(:, 1) - tr_mean) / tr_std;
G_all_std = single(cat(1, ...
    reshape(tr_std_all, 1, 1, []), ...
    reshape(sin(2 * ang_rad), 1, 1, []), ...
    reshape(cos(2 * ang_rad), 1, 1, [])));

%% Initialize

F = size(X_pad, 1);
nGlobal = size(G_all_std, 1);
params = init_hybrid_params(F, hiddenDim, cnnChannels, fusionDim, nComp, K, nGlobal);

fprintf('Hybrid params: mode=%s F=%d hidden=%d cnn=[%d %d %d] fusion=%d K=%d nComp=%d dropout=%.2f\n', ...
    modelMode, F, hiddenDim, cnnChannels(1), cnnChannels(2), cnnChannels(3), fusionDim, K, nComp, dropoutRate);

M = single(nodeMask);
if gpuUse
    fprintf('Moving params + loss weights to GPU. Bulk inputs transfer per batch.\n');
    params = dlupdate(@gpuArray, params);
    loss_w = gpuArray(loss_w);
    toGpu = @(x) gpuArray(x);
else
    toGpu = @(x) x;
end

%% Training

if training
    lossfcn = dlaccelerate(@hybrid_model_loss);
    clearCache(lossfcn);

    avgG = [];
    avgSqG = [];
    iteration = 0;
    bestValR2 = -inf;
    bestValLoss = inf;
    best_params = params;
    bestEpoch = 0;
    noImproveCount = 0;

    train_losses = zeros(maxEpochs, 1);
    val_losses = zeros(maxEpochs, 1);
    val_u3_R2 = nan(maxEpochs, 1);
    val_u3_mse = nan(maxEpochs, 1);

    numIters = max(1, ceil(nTrain / batchSize));
    tic
    for epoch = 1:maxEpochs
        lr = lr0 * decay_rate^(epoch - 1);
        perm = trainLocal(randperm(nTrain));
        epochLoss = 0;

        for b = 1:numIters
            si = (b - 1) * batchSize + 1;
            ei = min(b * batchSize, nTrain);
            idx = perm(si:ei);

            [loss, grad] = dlfeval(lossfcn, params, ...
                toGpu(X_pad(:, :, idx)), A_hat(idx), K, toGpu(M(:, :, idx)), ...
                toGpu(Dense_input(:, :, :, idx)), dlarray(toGpu(G_all_std(:, :, idx))), gpuUse, ...
                toGpu(Z_target(:, :, idx)), loss_w, dropoutRate, useGraph, useDense);

            grad = add_weight_decay_to_weights(grad, params, weight_decay);
            iteration = iteration + 1;
            [params, avgG, avgSqG] = adamupdate(params, grad, avgG, avgSqG, iteration, lr);

            epochLoss = epochLoss + double(gather(extractdata(loss)));
        end

        train_losses(epoch) = epochLoss / numIters;

        if epoch == 1 || mod(epoch, valFreq) == 0
            valLoss = evaluate_hybrid_loss(lossfcn, params, valLocal, evalBatchSize, ...
                X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, ...
                Z_target, loss_w, useGraph, useDense);
            val_losses(epoch) = valLoss;

            Zhat_val_std = predict_hybrid(params, valLocal, evalBatchSize, ...
                X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, useGraph, useDense);
            valMetrics = evaluate_hybrid_metrics(Zhat_val_std, Z_val, Z_mean, Z_std, ...
                coeff, u3_mean, U3_mat(val_idx, :), explained);
            val_u3_R2(epoch) = valMetrics.u3_R2;
            val_u3_mse(epoch) = valMetrics.u3_mse;

            if valMetrics.u3_R2 > bestValR2 + minDeltaR2
                bestValR2 = valMetrics.u3_R2;
                bestValLoss = valLoss;
                best_params = params;
                bestEpoch = epoch;
                noImproveCount = 0;
            else
                noImproveCount = noImproveCount + 1;
            end

            fprintf('Epoch %3d | train=%.5g | valLoss=%.5g | val_u3_R2=%.4f | lr=%.3e | patience=%d/%d\n', ...
                epoch, train_losses(epoch), valLoss, valMetrics.u3_R2, lr, noImproveCount, patience);

            if noImproveCount >= patience
                fprintf('Early stopping at epoch %d (best epoch=%d, val_u3_R2=%.4f)\n', ...
                    epoch, bestEpoch, bestValR2);
                break;
            end
        else
            fprintf('Epoch %3d | train=%.5g | lr=%.3e\n', epoch, train_losses(epoch), lr);
        end

        if mod(epoch, 2) == 0 || epoch == 1
            snap.epoch = epoch;
            snap.train_losses = train_losses;
            snap.val_losses = val_losses;
            snap.val_u3_R2 = val_u3_R2;
            snap.val_u3_mse = val_u3_mse;
            snap.params = dlupdate(@gather, params);
            save(progressFile, '-struct', 'snap', '-v7.3');
        end
    end
    toc

    params = best_params;
    cfg = struct('modelMode', modelMode, 'useGraph', useGraph, 'useDense', useDense, ...
        'K', K, 'hiddenDim', hiddenDim, 'cnnChannels', cnnChannels, 'fusionDim', fusionDim, ...
        'batchSize', batchSize, 'lr0', lr0, 'maxEpochs', maxEpochs, ...
        'bestEpoch', bestEpoch, 'bestValR2', bestValR2, 'bestValLoss', bestValLoss, ...
        'nComp', nComp, 'nGlobal', nGlobal);

    params_cpu = dlupdate(@gather, params);
    params_for_eval = params;
    params = params_cpu;
    save(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses', 'val_u3_R2', 'val_u3_mse');
    params = params_for_eval;
    fprintf('Saved best params (epoch %d, val_u3_R2=%.4f) -> %s\n', bestEpoch, bestValR2, paramFile);
else
    load(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses', 'val_u3_R2', 'val_u3_mse');
    if gpuUse, params = dlupdate(@gpuArray, params); end
    fprintf('Loaded params (best epoch=%d) from %s\n', cfg.bestEpoch, paramFile);
end

%% Final evaluation

Zhat_train_std = predict_hybrid(params, trainLocal, evalBatchSize, ...
    X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, useGraph, useDense);
Zhat_val_std = predict_hybrid(params, valLocal, evalBatchSize, ...
    X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, useGraph, useDense);
Zhat_test_std = predict_hybrid(params, testLocal, evalBatchSize, ...
    X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, useGraph, useDense);

trainMetrics = evaluate_hybrid_metrics(Zhat_train_std, Z_train, Z_mean, Z_std, ...
    coeff, u3_mean, U3_mat(train_idx, :), explained);
valMetrics = evaluate_hybrid_metrics(Zhat_val_std, Z_val, Z_mean, Z_std, ...
    coeff, u3_mean, U3_mat(val_idx, :), explained);
testMetrics = evaluate_hybrid_metrics(Zhat_test_std, Z_test, Z_mean, Z_std, ...
    coeff, u3_mean, U3_mat(test_idx, :), explained);

fprintf('\nHybrid vision final metrics (%s):\n', modelMode);
print_metrics('Train', trainMetrics, nComp);
print_metrics('Val', valMetrics, nComp);
print_metrics('Test', testMetrics, nComp);

metricsFile = fullfile(bestDir, 'final_metrics.mat');
save(metricsFile, 'trainMetrics', 'valMetrics', 'testMetrics', 'cfg');

%% Plots

valEpochs = find(val_losses > 0);
hLoss = figure;
plot(train_losses, '-r', 'LineWidth', 1.5); hold on;
plot(valEpochs, val_losses(valEpochs), '-b', 'LineWidth', 1.5);
xlabel('Epoch'); ylabel('Weighted MSE (standardized PCA)');
legend('Train', 'Validation', 'Location', 'northeast');
title(sprintf('Hybrid training curve (%s)', modelMode));
exportgraphics(hLoss, fullfile(bestDir, 'loss_curve.png'), 'Resolution', 150);

hR2 = figure;
plot(find(~isnan(val_u3_R2)), val_u3_R2(~isnan(val_u3_R2)), '-k', 'LineWidth', 1.5);
xlabel('Epoch'); ylabel('Validation u3 R2');
title(sprintf('Validation u3 R2 (%s)', modelMode));
exportgraphics(hR2, fullfile(bestDir, 'val_u3_r2.png'), 'Resolution', 150);

%% Plot: PCA component scatter (2 × 4)

Zhat_val = valMetrics.Zhat;   % nVal × nComp
hZ = figure('Position', [100 100 1200 500]);
for c = 1:nComp
    subplot(2, 4, c);
    scatter(Z_val(:,c), Zhat_val(:,c), 8, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    lo = min([Z_val(:,c); Zhat_val(:,c)]);
    hi = max([Z_val(:,c); Zhat_val(:,c)]);
    plot([lo hi], [lo hi], 'k--');
    xlabel('True'); ylabel('Pred');
    title(sprintf('PC%d  R2=%.3f', c, valMetrics.perPC_R2(c)));
    axis equal tight;
end
sgtitle(sprintf('PCA components: truth vs prediction (validation) [%s]', modelMode));
exportgraphics(hZ, fullfile(bestDir, 'pred_vs_truth.png'), 'Resolution', 150);

%% Plot: u3(s) profile overlays (3 samples)

U3_true_val = U3_mat(val_idx, :);
U3_pred_val = valMetrics.U3_pred;
nShow = min(3, nVal);
hU3 = figure('Position', [100 100 1100 320]);
for k = 1:nShow
    subplot(1, nShow, k);
    plot(s_grid, U3_true_val(k,:), 'k-',  'LineWidth', 1.5); hold on;
    plot(s_grid, U3_pred_val(k,:), 'r--', 'LineWidth', 1.5);
    xlabel('s'); ylabel('u_3 (m)');
    title(sprintf('Val sample %d', k));
    if k == 1, legend('True', 'Pred', 'Location', 'best'); end
end
sgtitle(sprintf('u_3(s) midpoint profile: truth vs prediction (validation) [%s]', modelMode));
exportgraphics(hU3, fullfile(bestDir, 'u3_profiles.png'), 'Resolution', 150);

%% Local helpers

function Zhat = predict_hybrid(params, localIdx, evalBatchSize, X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, useGraph, useDense)
    parts = cell(ceil(numel(localIdx) / evalBatchSize), 1);
    p = 0;
    for bi = 1:evalBatchSize:numel(localIdx)
        p = p + 1;
        idx = localIdx(bi:min(bi + evalBatchSize - 1, numel(localIdx)));
        parts{p} = hybrid_forward(params, toGpu(X_pad(:, :, idx)), A_hat(idx), K, ...
            toGpu(M(:, :, idx)), toGpu(Dense_input(:, :, :, idx)), ...
            dlarray(toGpu(G_all_std(:, :, idx))), gpuUse, 0.0, useGraph, useDense);
    end
    Zhat = cat(3, parts{:});
end

function avgLoss = evaluate_hybrid_loss(lossfcn, params, localIdx, evalBatchSize, X_pad, A_hat, K, M, Dense_input, G_all_std, gpuUse, toGpu, Z_target, loss_w, useGraph, useDense)
    totalLoss = 0;
    totalN = 0;
    for bi = 1:evalBatchSize:numel(localIdx)
        idx = localIdx(bi:min(bi + evalBatchSize - 1, numel(localIdx)));
        loss = dlfeval(lossfcn, params, ...
            toGpu(X_pad(:, :, idx)), A_hat(idx), K, toGpu(M(:, :, idx)), ...
            toGpu(Dense_input(:, :, :, idx)), dlarray(toGpu(G_all_std(:, :, idx))), gpuUse, ...
            toGpu(Z_target(:, :, idx)), loss_w, 0.0, useGraph, useDense);
        nBatch = numel(idx);
        totalLoss = totalLoss + double(gather(extractdata(loss))) * nBatch;
        totalN = totalN + nBatch;
    end
    avgLoss = totalLoss / max(totalN, 1);
end

function grad = add_weight_decay_to_weights(grad, params, wd)
    if wd == 0, return; end
    names = fieldnames(grad);
    for ii = 1:numel(names)
        name = names{ii};
        if isstruct(grad.(name))
            grad.(name) = add_weight_decay_to_weights(grad.(name), params.(name), wd);
        elseif isempty(grad.(name))
            grad.(name) = 0 * params.(name);
        elseif strcmp(name, 'W')
            grad.(name) = grad.(name) + wd * params.(name);
        end
    end
end

function print_metrics(label, metrics, nComp)
    fprintf('%s: u3_MSE=%.4g  u3_R2=%.4f  PCA_R2_w=%.4f  PCA_R2_unw=%.4f\n', ...
        label, metrics.u3_mse, metrics.u3_R2, metrics.R2_z_weighted, metrics.R2_z_unweighted);
    fprintf('  Per-PC R2:');
    for cc = 1:nComp
        fprintf(' PC%d=%.3f', cc, metrics.perPC_R2(cc));
    end
    fprintf('\n');
end
