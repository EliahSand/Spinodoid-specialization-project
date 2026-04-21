%% GNN to predict spinodoid sheet deformation (8-dim PCA target)

clear; close all; clc;
rng(1);

scriptDir = fileparts(mfilename('fullpath'));   % .../Matlab/GNN/models
gnnRoot   = fileparts(scriptDir);              % .../Matlab/GNN
addpath(genpath(fullfile(gnnRoot, 'helpers')));

%% Config

training     = true;
K            = 3;       % GCN layers (normally 3)
hiddenDim    = 64;      %64
poolDim      = 64;
batchSize    = 32;      %32
maxEpochs    = 120;       %120
lr0          = 3e-3;
decay_rate   = 0.99;
weight_decay = 1e-5;
valFreq      = 2;
patience     = 10;
minDelta     = 1e-4;
gpuUse       = canUseGPU;

runName    = datestr(now, 'yyyymmdd_HHMMSS');   % override for named runs, e.g. 'k3_h64'

targetsDir = fullfile(gnnRoot, 'data', 'dataset', 'targets');
samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');
bestDir    = fullfile(scriptDir, 'best', runName);
if ~isfolder(bestDir), mkdir(bestDir); end

paramFile    = fullfile(bestDir, 'best_params.mat');
progressFile = fullfile(bestDir, 'training_log.mat');

%% Load targets

t = load(fullfile(targetsDir, 'pca_targets.mat'), ...
    'Z_train', 'Z_val', 'Z_test', 'train_ids', 'val_ids', 'test_ids', 'n_components');
p = load(fullfile(targetsDir, 'pca_model.mat'),   'coeff', 'u3_mean', 's_grid');
u = load(fullfile(targetsDir, 'u3_targets.mat'),  'U3_mat');
sp = load(fullfile(targetsDir, 'split_indices.mat'), 'val_idx');

Z_train   = t.Z_train;    % n_train × nComp
Z_val     = t.Z_val;
Z_test    = t.Z_test;
train_ids = t.train_ids;
val_ids   = t.val_ids;
test_ids  = t.test_ids;
nComp     = t.n_components;

coeff   = p.coeff;     % N_GRID × nComp
u3_mean = p.u3_mean;   % 1 × N_GRID
s_grid  = p.s_grid;

U3_mat  = u.U3_mat;   % M × N_GRID
val_idx = sp.val_idx;

nTrain = numel(train_ids);
nVal   = numel(val_ids);
nTest  = numel(test_ids);

%% Load graphs

fprintf('Loading %d training graphs...\n',   nTrain);
[Xtr, eiTr, NTr] = load_graph_dataset(train_ids, samplesDir);
fprintf('Loading %d validation graphs...\n', nVal);
[Xvl, eiVl, NVl] = load_graph_dataset(val_ids, samplesDir);
fprintf('Loading %d test graphs...\n',       nTest);
[Xts, eiTs, NTs] = load_graph_dataset(test_ids, samplesDir);

% Stack all splits; track local index ranges
allX  = [Xtr;  Xvl;  Xts];
allEI = [eiTr; eiVl; eiTs];
allN  = [NTr;  NVl;  NTs];
nAll  = nTrain + nVal + nTest;

trainLocal = (1 : nTrain)';
valLocal   = (nTrain+1       : nTrain+nVal)';
testLocal  = (nTrain+nVal+1  : nAll)';

%% Normalize features & build adjacency

trainMask = false(nAll, 1);
trainMask(trainLocal) = true;

maxN = max(allN);
fprintf('maxN=%d  nAll=%d\n', maxN, nAll);

fprintf('Normalizing node features...\n');
[X_pad, nodeMask, normStats] = pad_and_normalize_graphs(allX, allN, trainMask);
% X_pad:    3 × maxN × nAll  (single)
% nodeMask: 1 × maxN × nAll  (logical)

fprintf('Building normalized adjacency matrices...\n');
A_hat = build_norm_adjacency(allEI, allN, maxN);

%% Standardize PCA targets (per-component, train stats only)

Z_mean = mean(Z_train, 1);             % 1 × nComp
Z_std  = std(Z_train, 0, 1);
Z_std(Z_std < 1e-8) = 1e-8;

normStats.Z_mean = Z_mean;
normStats.Z_std  = Z_std;

Z_all     = [Z_train; Z_val; Z_test];                  % nAll × nComp
Z_all_std = (Z_all - Z_mean) ./ Z_std;

% Target tensor: nComp × 1 × nAll
Z_target = single(reshape(Z_all_std', nComp, 1, []));

%% Initialize model

F      = 3;   % node features: x, y, radius
params = init_params(F, hiddenDim, poolDim, nComp, K);

fprintf('Params: F=%d  hidden=%d  pool=%d  K=%d  nComp=%d\n', F, hiddenDim, poolDim, K, nComp);

%% Training

M = single(nodeMask);   % 1 × maxN × nAll  (used as mask in loss)

if training
    lossfcn = dlaccelerate(@model_loss);
    clearCache(lossfcn);

    avgG = []; avgSqG = []; iteration = 0;
    bestValLoss = inf; best_params = params; bestEpoch = 0; noImproveCount = 0;

    train_losses = zeros(maxEpochs, 1);
    val_losses   = zeros(maxEpochs, 1);

    numIters = max(1, ceil(nTrain / batchSize));

    tic
    for epoch = 1:maxEpochs
        lr   = lr0 * decay_rate^(epoch - 1);
        perm = trainLocal(randperm(nTrain));
        epochLoss = 0;

        for b = 1:numIters
            si  = (b-1)*batchSize + 1;
            ei  = min(b*batchSize, nTrain);
            idx = perm(si:ei);

            [loss, grad] = dlfeval(lossfcn, params, ...
                X_pad(:,:,idx), A_hat(idx), K, M(:,:,idx), gpuUse, Z_target(:,:,idx));

            grad      = dlupdate(@(g,p) g + weight_decay*p, grad, params);
            iteration = iteration + 1;
            [params, avgG, avgSqG] = adamupdate(params, grad, avgG, avgSqG, iteration, lr);

            epochLoss = epochLoss + double(gather(extractdata(loss)));
        end

        train_losses(epoch) = epochLoss / numIters;

        if epoch == 1 || mod(epoch, valFreq) == 0
            [valLoss, ~] = dlfeval(lossfcn, params, ...
                X_pad(:,:,valLocal), A_hat(valLocal), K, M(:,:,valLocal), gpuUse, Z_target(:,:,valLocal));
            valLoss = double(gather(extractdata(valLoss)));
            val_losses(epoch) = valLoss;

            if valLoss < bestValLoss - minDelta
                bestValLoss = valLoss;
                best_params = params;
                bestEpoch   = epoch;
                noImproveCount = 0;
            else
                noImproveCount = noImproveCount + 1;
            end

            fprintf('Epoch %3d | train=%.5g | val=%.5g | lr=%.3e | patience=%d/%d\n', ...
                epoch, train_losses(epoch), valLoss, lr, noImproveCount, patience);

            if noImproveCount >= patience
                fprintf('Early stopping at epoch %d (best epoch=%d, val=%.5g)\n', epoch, bestEpoch, bestValLoss);
                break;
            end
        else
            fprintf('Epoch %3d | train=%.5g | lr=%.3e\n', epoch, train_losses(epoch), lr);
        end

        if mod(epoch, 2) == 0 || epoch == 1
            save(progressFile, 'epoch', 'train_losses', 'val_losses', 'params', '-v7.3');
        end
    end
    toc

    params = best_params;
    cfg = struct('K', K, 'hiddenDim', hiddenDim, 'poolDim', poolDim, ...
                 'batchSize', batchSize, 'lr0', lr0, 'maxEpochs', maxEpochs, ...
                 'bestEpoch', bestEpoch, 'nComp', nComp);
    save(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses');
    fprintf('Saved best params (epoch %d) → %s\n', bestEpoch, paramFile);

else
    if ~isfile(paramFile)
        error('Set training=true or provide %s', paramFile);
    end
    load(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses');
    fprintf('Loaded params (best epoch=%d) from %s\n', cfg.bestEpoch, paramFile);
end

%% Evaluate on validation set

Zhat_std = gcn_forward(params, X_pad(:,:,valLocal), A_hat(valLocal), K, M(:,:,valLocal), gpuUse);
Zhat_std = gather(extractdata(Zhat_std));   % nComp × 1 × nVal
Zhat_val = squeeze(Zhat_std).' .* Z_std + Z_mean;   % nVal × nComp

% Z-space metrics
err_z  = Zhat_val - Z_val;
mse_z  = mean(err_z(:).^2);
R2_z   = 1 - sum(err_z(:).^2) / sum((Z_val(:) - mean(Z_val(:))).^2);
fprintf('\nValidation Z-space:   MSE=%.4g  R2=%.4f\n', mse_z, R2_z);

% u3-space metrics (decode PCA)
U3_pred = Zhat_val * coeff' + u3_mean;      % nVal × N_GRID
U3_true = U3_mat(val_idx, :);

err_u3 = U3_pred - U3_true;
mse_u3 = mean(err_u3(:).^2);
R2_u3  = 1 - sum(err_u3(:).^2) / sum((U3_true(:) - mean(U3_true(:))).^2);
fprintf('Validation u3-space:  MSE=%.4g  R2=%.4f\n\n', mse_u3, R2_u3);

%% Plot: loss curve

valEpochs = find(val_losses > 0);
hLoss = figure;
plot(train_losses, '-r', 'LineWidth', 1.5); hold on;
plot(valEpochs, val_losses(valEpochs), '-b', 'LineWidth', 1.5);
xlabel('Epoch'); ylabel('MSE (standardized Z)');
legend('Train', 'Validation', 'Location', 'northeast');
title('Training curve');
exportgraphics(hLoss, fullfile(bestDir, 'loss_curve.png'), 'Resolution', 150);

%% Plot: PCA component scatter (2 × 4)

hZ = figure('Position', [100 100 1200 500]);
for c = 1:nComp
    subplot(2, 4, c);
    scatter(Z_val(:,c), Zhat_val(:,c), 8, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    lo = min([Z_val(:,c); Zhat_val(:,c)]);
    hi = max([Z_val(:,c); Zhat_val(:,c)]);
    plot([lo hi], [lo hi], 'k--');
    xlabel('True'); ylabel('Pred');
    title(sprintf('PC%d', c));
    axis equal tight;
end
sgtitle('PCA components: truth vs prediction (validation)');
exportgraphics(hZ, fullfile(bestDir, 'pred_vs_truth.png'), 'Resolution', 150);

%% Plot: u3(s) profile overlays (3 samples)

nShow = min(3, nVal);
hU3 = figure('Position', [100 100 1100 320]);
for k = 1:nShow
    subplot(1, nShow, k);
    plot(s_grid, U3_true(k,:), 'k-',  'LineWidth', 1.5); hold on;
    plot(s_grid, U3_pred(k,:), 'r--', 'LineWidth', 1.5);
    xlabel('s'); ylabel('u_3 (m)');
    title(sprintf('Val sample %d', k));
    if k == 1, legend('True', 'Pred', 'Location', 'best'); end
end
sgtitle('u_3(s) midpoint profile: truth vs prediction');
exportgraphics(hU3, fullfile(bestDir, 'u3_profiles.png'), 'Resolution', 150);
