% step7_fit_pca_u3.m
%
% Fits PCA on the TRAINING split of u3(s) profiles and projects all splits.
% PCA is fit on train only to prevent data leakage into val/test.
%
% Inputs:   data/dataset/targets/u3_targets.mat
%           data/dataset/targets/split_indices.mat
% Outputs:  data/dataset/targets/pca_model.mat      — PCA basis (train only)
%           data/dataset/targets/pca_targets.mat    — projections for all splits
%           data/dataset/targets/pca_manifest.json

%% ---- Config ---------------------------------------------------------------
REPO_ROOT   = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'targets');

N_COMPONENTS = 8;    % 5 PCA components gives 99% variance, fewest GNN parameters to predict for the GNN
VAR_THRESHOLD = 0.99; % logged only — not used to truncate

%% ---- Load -----------------------------------------------------------------
targetsFile = fullfile(TARGETS_DIR, 'u3_targets.mat');
splitFile   = fullfile(TARGETS_DIR, 'split_indices.mat');

if ~isfile(targetsFile), error('step7:missing', 'Run step5 first.'); end
if ~isfile(splitFile),   error('step7:missing', 'Run step6 first.'); end

fprintf('Loading targets and split...\n');
t = load(targetsFile, 'U3_mat', 's_grid', 'sample_ids');
s = load(splitFile,   'train_idx', 'val_idx', 'test_idx', 'rng_seed');

U3        = t.U3_mat;       % M x N_GRID
M         = size(U3, 1);
train_idx = s.train_idx;
val_idx   = s.val_idx;
test_idx  = s.test_idx;

fprintf('M=%d  train=%d  val=%d  test=%d\n', M, numel(train_idx), numel(val_idx), numel(test_idx));

%% ---- Fit PCA on train only ------------------------------------------------
fprintf('Fitting PCA on %d training samples (n_components=%d)...\n', numel(train_idx), N_COMPONENTS);

U3_train = U3(train_idx, :);
u3_mean  = mean(U3_train, 1);          % 1 x N_GRID
U3_centered = U3_train - u3_mean;

[coeff, ~, ~, ~, explained] = pca(U3_centered, 'NumComponents', N_COMPONENTS);
% coeff: N_GRID x N_COMPONENTS

cum_var = cumsum(explained);
n_for_threshold = find(cum_var >= VAR_THRESHOLD * 100, 1);
fprintf('Variance explained: ');
fprintf('PC%d=%.2f%%  ', [1:min(5,N_COMPONENTS); explained(1:min(5,N_COMPONENTS))']);
fprintf('\n');
fprintf('Cumulative %.0f%% variance reached at %d components.\n', VAR_THRESHOLD*100, n_for_threshold);

%% ---- Project all splits ---------------------------------------------------
project = @(X) (X - u3_mean) * coeff;   % (n x N_GRID) -> (n x N_COMPONENTS)

Z_train = project(U3(train_idx, :));
Z_val   = project(U3(val_idx,   :));
Z_test  = project(U3(test_idx,  :));

%% ---- Reconstruction error (sanity check) ----------------------------------
U3_recon_train = Z_train * coeff' + u3_mean;
recon_err_train = mean(mean((U3(train_idx,:) - U3_recon_train).^2, 2));
fprintf('Mean squared reconstruction error (train): %.3e m^2\n', recon_err_train);

%% ---- Save PCA model -------------------------------------------------------
created_at     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
schema_version = 1;
n_components   = N_COMPONENTS;
s_grid         = t.s_grid;

pcaModelFile = fullfile(TARGETS_DIR, 'pca_model.mat');
save(pcaModelFile, 'coeff', 'u3_mean', 'explained', 'n_components', ...
    's_grid', 'schema_version', 'created_at');
fprintf('Saved PCA model: %s\n', pcaModelFile);

%% ---- Save projected targets -----------------------------------------------
sample_ids = t.sample_ids;
train_ids  = sample_ids(train_idx);
val_ids    = sample_ids(val_idx);
test_ids   = sample_ids(test_idx);

pcaTargetsFile = fullfile(TARGETS_DIR, 'pca_targets.mat');
save(pcaTargetsFile, ...
    'Z_train', 'Z_val', 'Z_test', ...
    'train_ids', 'val_ids', 'test_ids', ...
    'train_idx', 'val_idx', 'test_idx', ...
    'n_components', 'schema_version', 'created_at');
fprintf('Saved PCA targets: %s\n', pcaTargetsFile);

%% ---- Manifest -------------------------------------------------------------
try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

manifest = struct( ...
    'n_components',           N_COMPONENTS, ...
    'var_threshold_pct',      VAR_THRESHOLD * 100, ...
    'n_components_for_threshold', n_for_threshold, ...
    'explained_pct',          explained', ...
    'recon_mse_train',        recon_err_train, ...
    'pca_fit_on',             'train_only', ...
    'n_train',                numel(train_idx), ...
    'n_val',                  numel(val_idx), ...
    'n_test',                 numel(test_idx), ...
    'split_rng_seed',         s.rng_seed, ...
    'schema_version',         schema_version, ...
    'git_sha',                git_sha, ...
    'created_at',             created_at);

manifestFile = fullfile(TARGETS_DIR, 'pca_manifest.json');
fid = fopen(manifestFile, 'w');
fprintf(fid, '%s\n', jsonencode(manifest));
fclose(fid);
fprintf('Manifest: %s\n', manifestFile);
fprintf('Done.\n');
