%% GNN to predict spinodoid sheet deformation (8-dim PCA target)

clear; close all; clc;
rng(1);

scriptDir = fileparts(mfilename('fullpath'));   % .../Matlab/GNN/models
gnnRoot   = fileparts(scriptDir);              % .../Matlab/GNN
addpath(genpath(fullfile(gnnRoot, 'helpers')));

%% Config

training       = true;
subsetFraction = 0.5;     % 1.0 = full dataset; 0.1 = 10% feks
K              = 12;       % GCN layers (normally 3)
hiddenDim      = 128;     % GCN hidden dimension
poolDim        = 128;     % (unused, replaced by readoutDim)
readoutDim     = 128;     % Readout head dimension after attention pooling
batchSize      = 32;      %32
maxEpochs      = 120;       %120
lr0            = 1e-3;
decay_rate     = 0.99;
weight_decay   = 1e-5;
valFreq        = 5;
patience       = 10;
minDelta       = 1e-4;
gpuUse         = canUseGPU;

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

%% Optional subset (fast-iteration mode)

if subsetFraction < 1.0
    assert(subsetFraction > 0, 'subsetFraction must be > 0');
    nTrSub = max(1, round(numel(train_ids) * subsetFraction));
    nVSub  = max(1, round(numel(val_ids)   * subsetFraction));
    nTsSub = max(1, round(numel(test_ids)  * subsetFraction));

    trSel = sort(randperm(numel(train_ids), nTrSub));
    vSel  = sort(randperm(numel(val_ids),   nVSub));
    tsSel = sort(randperm(numel(test_ids),  nTsSub));

    train_ids = train_ids(trSel);  Z_train = Z_train(trSel, :);
    val_ids   = val_ids(vSel);     Z_val   = Z_val(vSel, :);   val_idx = val_idx(vSel);
    test_ids  = test_ids(tsSel);   Z_test  = Z_test(tsSel, :);

    fprintf('SUBSET: using %.0f%% -> train=%d val=%d test=%d\n', ...
        100*subsetFraction, nTrSub, nVSub, nTsSub);
end

nTrain = numel(train_ids);
nVal   = numel(val_ids);
nTest  = numel(test_ids);

%% Load graphs

fprintf('Loading %d training graphs...\n',   nTrain);
[Xtr, eiTr, NTr, Gtr] = load_graph_dataset(train_ids, samplesDir);
fprintf('Loading %d validation graphs...\n', nVal);
[Xvl, eiVl, NVl, Gvl] = load_graph_dataset(val_ids, samplesDir);
fprintf('Loading %d test graphs...\n',       nTest);
[Xts, eiTs, NTs, Gts] = load_graph_dataset(test_ids, samplesDir);

% Stack all splits; track local index ranges
allX  = [Xtr;  Xvl;  Xts];
allEI = [eiTr; eiVl; eiTs];
allN  = [NTr;  NVl;  NTs];
G_all = [Gtr;  Gvl;  Gts];   % nAll × 2  ([tr_ratio, ang_deg])
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
[X_pad, nodeMask, normStats] = pad_and_normalize_graphs(allX, allN, trainMask, []);
% X_pad:    4 × maxN × nAll  (single)  [x, y, radius, boundary]
% nodeMask: 1 × maxN × nAll  (logical)

fprintf('Building normalized adjacency matrices...\n');
A_hat = build_norm_adjacency(allEI, allN, maxN);

%% Standardize PCA targets (per-component, train stats only)

Z_mean = mean(Z_train, 1);             % 1 × nComp
Z_std  = std(Z_train, 0, 1);
Z_std(Z_std < 1e-8) = 1e-8;

normStats.Z_mean = Z_mean;
normStats.Z_std  = Z_std;

% Variance-proportional loss weights: w_i ∝ Var(PC_i)
% Weighted MSE on standardized Z = raw-Z MSE = u3-space MSE (orthonormal basis)
loss_w = single(Z_std.^2);
loss_w = loss_w / mean(loss_w);          % normalize: mean weight = 1
loss_w = reshape(loss_w, nComp, 1, 1);  % nComp × 1 × 1
normStats.loss_w = loss_w;

Z_all     = [Z_train; Z_val; Z_test];                  % nAll × nComp
Z_all_std = (Z_all - Z_mean) ./ Z_std;

% Target tensor: nComp × 1 × nAll
Z_target = single(reshape(Z_all_std', nComp, 1, []));

%% Standardize global graph features
G_mean = mean(G_all(trainLocal, :), 1);        % 1 × 2
G_std  = std(G_all(trainLocal, :), 0, 1);
G_std(G_std < 1e-8) = 1e-8;

normStats.G_mean = G_mean;
normStats.G_std  = G_std;

G_all_std = single(reshape(((G_all - G_mean) ./ G_std)', 2, 1, []));  % 2 × 1 × nAll

%% Initialize model

F          = 4;     % node features: x, y, radius, boundary
nGlobal    = 2;     % global features: tr_ratio, ang_deg (concatenated after pooling)
params     = init_params(F, hiddenDim, poolDim, nComp, K, nGlobal, readoutDim);

fprintf('Params: F=%d  hidden=%d  readout=%d  K=%d  nComp=%d  nGlobal=%d\n', ...
    F, hiddenDim, readoutDim, K, nComp, nGlobal);

%% Move data and parameters to GPU (when available)

M = single(nodeMask);   % 1 × maxN × nAll  (used as mask in loss)

if gpuUse
    fprintf('Moving tensors to GPU...\n');
    X_pad     = gpuArray(X_pad);
    M         = gpuArray(M);
    Z_target  = gpuArray(Z_target);
    G_all_std = gpuArray(G_all_std);
    loss_w    = gpuArray(loss_w);
    params    = dlupdate(@gpuArray, params);
end

%% Training

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
                X_pad(:,:,idx), A_hat(idx), K, M(:,:,idx), gpuUse, Z_target(:,:,idx), loss_w, ...
                dlarray(G_all_std(:,:,idx)));

            grad      = dlupdate(@(g,p) g + weight_decay*p, grad, params);
            iteration = iteration + 1;
            [params, avgG, avgSqG] = adamupdate(params, grad, avgG, avgSqG, iteration, lr);

            epochLoss = epochLoss + double(gather(extractdata(loss)));
        end

        train_losses(epoch) = epochLoss / numIters;

        if epoch == 1 || mod(epoch, valFreq) == 0
            [valLoss, ~] = dlfeval(lossfcn, params, ...
                X_pad(:,:,valLocal), A_hat(valLocal), K, M(:,:,valLocal), gpuUse, ...
                Z_target(:,:,valLocal), loss_w, dlarray(G_all_std(:,:,valLocal)));
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
            snap.epoch        = epoch;
            snap.train_losses = train_losses;
            snap.val_losses   = val_losses;
            snap.params       = dlupdate(@gather, params);
            save(progressFile, '-struct', 'snap', '-v7.3');
        end
    end
    toc

    params = best_params;
    cfg = struct('K', K, 'hiddenDim', hiddenDim, 'readoutDim', readoutDim, ...
                 'batchSize', batchSize, 'lr0', lr0, 'maxEpochs', maxEpochs, ...
                 'bestEpoch', bestEpoch, 'nComp', nComp, 'nGlobal', nGlobal);

    % Save on CPU for portability; keep GPU copy for downstream evaluation
    params_cpu = dlupdate(@gather, params);
    params_for_save = params;
    params = params_cpu; 
    save(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses');
    params = params_for_save;
    fprintf('Saved best params (epoch %d) → %s\n', bestEpoch, paramFile);

else
    if ~isfile(paramFile)
        error('Set training=true or provide %s', paramFile);
    end
    load(paramFile, 'params', 'cfg', 'normStats', 'train_losses', 'val_losses');
    if gpuUse
        params = dlupdate(@gpuArray, params);
    end
    fprintf('Loaded params (best epoch=%d) from %s\n', cfg.bestEpoch, paramFile);
end

%% Evaluate on validation set

Zhat_std = gcn_forward(params, X_pad(:,:,valLocal), A_hat(valLocal), K, M(:,:,valLocal), gpuUse, ...
    dlarray(G_all_std(:,:,valLocal)));
Zhat_std = gather(extractdata(Zhat_std));   % nComp × 1 × nVal
Zhat_val = squeeze(Zhat_std).' .* Z_std + Z_mean;   % nVal × nComp

% Z-space metrics
err_z    = Zhat_val - Z_val;                               % nVal × nComp
ss_res_c = sum(err_z.^2, 1);                              % 1 × nComp
ss_tot_c = sum((Z_val - mean(Z_val, 1)).^2, 1);

% Per-component R²
perPC_R2 = 1 - ss_res_c ./ max(ss_tot_c, eps);

% Unweighted aggregate (treats all PCs equally — for comparison with old runs)
R2_z_unw = 1 - sum(ss_res_c) / sum(ss_tot_c);

% Variance-weighted aggregate (w_i ∝ Var(PC_i) → equivalent to u3-space R²)
w_var      = Z_std.^2;
R2_z_wt    = 1 - sum(w_var .* ss_res_c) / sum(w_var .* ss_tot_c);

fprintf('\nValidation Z-space:\n');
fprintf('  Unweighted R²:        %.4f   (equal weight per PC)\n', R2_z_unw);
fprintf('  Variance-weighted R²: %.4f   (proportional to PC variance, ≈ u3-space R²)\n', R2_z_wt);
fprintf('  Per-PC R²:');
for c = 1:nComp
    fprintf('  PC%d=%.3f', c, perPC_R2(c));
end
fprintf('\n');

% u3-space metrics (decode PCA)
U3_pred = Zhat_val * coeff' + u3_mean;      % nVal × N_GRID
U3_true = U3_mat(val_idx, :);

err_u3 = U3_pred - U3_true;
mse_u3 = mean(err_u3(:).^2);
R2_u3  = 1 - sum(err_u3(:).^2) / sum((U3_true(:) - mean(U3_true(:))).^2);
fprintf('Validation u3-space:  MSE=%.4g  R²=%.4f\n\n', mse_u3, R2_u3);

%% Plot: loss curve

valEpochs = find(val_losses > 0);
hLoss = figure;
plot(train_losses, '-r', 'LineWidth', 1.5); hold on;
plot(valEpochs, val_losses(valEpochs), '-b', 'LineWidth', 1.5);
xlabel('Epoch'); ylabel('Variance-weighted MSE (standardized Z)');
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
