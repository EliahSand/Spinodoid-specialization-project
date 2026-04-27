% step8_verify_pca_target.m
%
% Diagnostics to verify that PCA coefficients are a suitable GNN target.
% Answers six questions; prints a PASS/WARN verdict for each and saves
% results to data/dataset/targets/diagnostics/.
%
% Run order: step6 -> step7 -> step8
%
% Inputs:   u3_targets.mat, split_indices.mat, pca_model.mat, pca_targets.mat
% Outputs:  diagnostics/pca_diagnostics.json
%           diagnostics/scree.png
%           diagnostics/coeff_hist.png
%           diagnostics/worst_recon/
%           diagnostics/best_recon/
%           diagnostics/avg_recon/

%% ---- Config ---------------------------------------------------------------
REPO_ROOT   = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'targets');
DIAG_DIR    = fullfile(TARGETS_DIR, 'diagnostics');

N_WORST      = 5;     % samples to plot for worst-reconstruction panel
N_BEST       = 5;     % samples to plot for best-reconstruction panel
N_AVG        = 5;     % samples to plot for average-reconstruction panel
N_HIST_PCS   = 8;     % PCs to show in coefficient histogram
STABILITY_SEED = 43;  % alternate seed for basis stability check (Q3)

% Pass thresholds (logged in JSON; not hard errors)
THRESH_RELATIVE_RMSE_RATIO = 2.0;   % val/test vs train relative RMSE
THRESH_K_99PCT             = 32;    % components for 99% variance
THRESH_SUBSPACE_DEG        = 10.0;  % max subspace angle (degrees)

if ~isfolder(DIAG_DIR), mkdir(DIAG_DIR); end

%% ---- Load -----------------------------------------------------------------
fprintf('Loading data...\n');
t  = load(fullfile(TARGETS_DIR, 'u3_targets.mat'));
sp = load(fullfile(TARGETS_DIR, 'split_indices.mat'));
pm = load(fullfile(TARGETS_DIR, 'pca_model.mat'));
pt = load(fullfile(TARGETS_DIR, 'pca_targets.mat'));

U3         = t.U3_mat;          % M x N_GRID
s_grid     = t.s_grid;
sample_ids = t.sample_ids;
coeff      = pm.coeff;           % N_GRID x K
u3_mean    = pm.u3_mean;         % 1 x N_GRID
explained  = pm.explained;       % K x 1
K          = size(coeff, 2);

train_idx = sp.train_idx;
val_idx   = sp.val_idx;
test_idx  = sp.test_idx;
Z_train   = pt.Z_train;          % n_train x K
Z_val     = pt.Z_val;
Z_test    = pt.Z_test;

fprintf('K=%d components  train=%d  val=%d  test=%d\n', ...
    K, numel(train_idx), numel(val_idx), numel(test_idx));

diag = struct();  % collects all numeric results

%% =========================================================================
%% Q1 — Reconstruction fidelity (basis generalizes to held-out data?)
%% =========================================================================
fprintf('\n--- Q1: Reconstruction fidelity ---\n');

function [abs_rmse, rel_rmse, worst_rel, worst_id] = recon_metrics(Z, U3_true, coeff, u3_mean, sample_ids)
    U3_recon = Z * coeff' + u3_mean;
    err      = U3_true - U3_recon;
    per_sample_rmse = sqrt(mean(err.^2, 2));
    per_sample_max  = max(abs(U3_true), [], 2);
    per_sample_rel  = per_sample_rmse ./ per_sample_max;
    abs_rmse  = mean(per_sample_rmse);
    rel_rmse  = mean(per_sample_rel);
    [worst_rel, wi] = max(per_sample_rel);
    worst_id  = sample_ids{wi};
end

[q1_train_abs, q1_train_rel, ~, ~]              = recon_metrics(Z_train, U3(train_idx,:), coeff, u3_mean, sample_ids(train_idx));
[q1_val_abs,   q1_val_rel,   q1_val_worst_rel,   q1_val_worst_id]   = recon_metrics(Z_val,   U3(val_idx,:),   coeff, u3_mean, sample_ids(val_idx));
[q1_test_abs,  q1_test_rel,  q1_test_worst_rel,  q1_test_worst_id]  = recon_metrics(Z_test,  U3(test_idx,:),  coeff, u3_mean, sample_ids(test_idx));

fprintf('  Abs RMSE  — train: %.3e m   val: %.3e m   test: %.3e m\n', q1_train_abs, q1_val_abs, q1_test_abs);
fprintf('  Rel RMSE  — train: %.3f%%   val: %.3f%%   test: %.3f%%\n', 100*q1_train_rel, 100*q1_val_rel, 100*q1_test_rel);
fprintf('  Worst val: %s  (%.3f%%)\n', q1_val_worst_id, 100*q1_val_worst_rel);

q1_pass = (q1_val_rel / q1_train_rel) < THRESH_RELATIVE_RMSE_RATIO && ...
          (q1_test_rel / q1_train_rel) < THRESH_RELATIVE_RMSE_RATIO;
fprintf('  Q1: %s\n', verdict(q1_pass));

diag.q1_train_abs_rmse   = q1_train_abs;
diag.q1_train_rel_rmse   = q1_train_rel;
diag.q1_val_abs_rmse     = q1_val_abs;
diag.q1_val_rel_rmse     = q1_val_rel;
diag.q1_test_abs_rmse    = q1_test_abs;
diag.q1_test_rel_rmse    = q1_test_rel;
diag.q1_val_worst_rel    = q1_val_worst_rel;
diag.q1_val_worst_id     = q1_val_worst_id;
diag.q1_test_worst_rel   = q1_test_worst_rel;
diag.q1_test_worst_id    = q1_test_worst_id;
diag.q1_pass             = q1_pass;

%% =========================================================================
%% Q2 — Scree / variance curve
%% =========================================================================
fprintf('\n--- Q2: Scree / variance ---\n');

cum_var = cumsum(explained(1:K));
thresholds = [90 95 99 99.9];
k_at = zeros(size(thresholds));
for ti = 1:numel(thresholds)
    idx = find(cum_var >= thresholds(ti), 1);
    if isempty(idx), idx = K; end
    k_at(ti) = idx;
    fprintf('  k for %.1f%% variance: %d\n', thresholds(ti), idx);
end
k_99 = k_at(3);
q2_pass = k_99 <= THRESH_K_99PCT;
fprintf('  Q2: %s\n', verdict(q2_pass));

diag.q2_k_for_90pct  = k_at(1);
diag.q2_k_for_95pct  = k_at(2);
diag.q2_k_for_99pct  = k_at(3);
diag.q2_k_for_999pct = k_at(4);
diag.q2_pass         = q2_pass;

% Plot
fig = figure('Visible','off');
plot(1:K, cum_var, 'b-o', 'MarkerSize', 4); hold on;
yline(99, 'r--', '99%');
yline(95, 'k--', '95%');
xlabel('Number of components'); ylabel('Cumulative variance explained (%)');
title('PCA scree — u3(s) profiles (train fit)');
grid on;
exportgraphics(fig, fullfile(DIAG_DIR, 'scree.png'), 'Resolution', 150);
close(fig);
fprintf('  Saved scree.png\n');

%% =========================================================================
%% Q3 — Basis stability (refit on alternate seed)
%% =========================================================================
fprintf('\n--- Q3: Basis stability ---\n');

% Rebuild alternate train split (same stratification, different seed)
rng(STABILITY_SEED);
tr_vals  = unique(t.tr_ratio);
ang_vals = unique(t.ang_deg);
alt_train_idx = [];
for ti2 = 1:numel(tr_vals)
    for ai2 = 1:numel(ang_vals)
        bin_mask = t.tr_ratio == tr_vals(ti2) & t.ang_deg == ang_vals(ai2);
        bin_idx  = find(bin_mask);
        if isempty(bin_idx), continue; end
        bin_idx = bin_idx(randperm(numel(bin_idx)));
        n = numel(bin_idx);
        n_train2 = max(1, round(n * sp.train_frac));
        alt_train_idx = [alt_train_idx; bin_idx(1:n_train2)]; %#ok<AGROW>
    end
end

U3_alt = U3(alt_train_idx, :);
u3_mean_alt = mean(U3_alt, 1);
[coeff_alt, ~, ~, ~, ~] = pca(U3_alt - u3_mean_alt, 'NumComponents', K);

% Subspace angle between top-K subspaces (radians -> degrees)
ang_rad = subspace(coeff, coeff_alt);
ang_deg_val = rad2deg(ang_rad);
fprintf('  Subspace angle between seed-42 and seed-43 bases: %.2f deg\n', ang_deg_val);

% Per-component cosine similarity (align signs first)
cos_sim = zeros(K,1);
for k = 1:K
    c1 = coeff(:,k); c2 = coeff_alt(:,k);
    c2 = c2 * sign(c1'*c2);   % align sign
    cos_sim(k) = abs(c1'*c2) / (norm(c1)*norm(c2));
end
fprintf('  Mean per-component cosine similarity: %.4f\n', mean(cos_sim));

q3_pass = ang_deg_val < THRESH_SUBSPACE_DEG;
fprintf('  Q3: %s\n', verdict(q3_pass));

diag.q3_subspace_angle_deg      = ang_deg_val;
diag.q3_mean_cosine_similarity  = mean(cos_sim);
diag.q3_pass                    = q3_pass;

%% =========================================================================
%% Q4 — Coefficient distributions
%% =========================================================================
fprintf('\n--- Q4: Coefficient distributions ---\n');

n_show = min(N_HIST_PCS, K);
fig = figure('Visible','off');
fig.Position = [0 0 1200 600];
for k = 1:n_show
    subplot(2, ceil(n_show/2), k);
    histogram(Z_train(:,k), 40);
    title(sprintf('PC%d  \\sigma=%.3f', k, std(Z_train(:,k))));
    xlabel('Coefficient'); ylabel('Count');
end
sgtitle('PCA coefficient distributions (train)');
exportgraphics(fig, fullfile(DIAG_DIR, 'coeff_hist.png'), 'Resolution', 150);
close(fig);

coeff_stats = struct();
for k = 1:K
    coeff_stats.(sprintf('pc%d_mean', k))     = mean(Z_train(:,k));
    coeff_stats.(sprintf('pc%d_std', k))      = std(Z_train(:,k));
    coeff_stats.(sprintf('pc%d_skewness', k)) = skewness(Z_train(:,k));
end
diag.q4_coeff_stats = coeff_stats;
fprintf('  Saved coeff_hist.png\n');
fprintf('  PC1 std=%.4f  PC2 std=%.4f  PC3 std=%.4f\n', ...
    std(Z_train(:,1)), std(Z_train(:,min(2,K))), std(Z_train(:,min(3,K))));

%% =========================================================================
%% Q5 — Worst / average / best reconstructions (visual)
%% =========================================================================
fprintf('\n--- Q5: Worst / average / best reconstructions ---\n');

U3_val_recon = Z_val * coeff' + u3_mean;
err_val = U3(val_idx,:) - U3_val_recon;
per_sample_rel_val = sqrt(mean(err_val.^2, 2)) ./ max(abs(U3(val_idx,:)), [], 2);
[~, sort_ord] = sort(per_sample_rel_val, 'descend');   % sort_ord(1) = worst
n_val = numel(sort_ord);

worst_k = sort_ord(1 : min(N_WORST, n_val));
best_k  = sort_ord(n_val : -1 : max(1, n_val - N_BEST + 1));   % rank 1 = best
mid     = round(n_val / 2);
half    = floor(N_AVG / 2);
avg_range = max(1, mid - half) : min(n_val, mid - half + N_AVG - 1);
avg_k   = sort_ord(avg_range);

plot_recon_batch(fullfile(DIAG_DIR, 'worst_recon'), 'Worst', worst_k, per_sample_rel_val, val_idx, U3, U3_val_recon, s_grid, sample_ids);
plot_recon_batch(fullfile(DIAG_DIR, 'best_recon'),  'Best',  best_k,  per_sample_rel_val, val_idx, U3, U3_val_recon, s_grid, sample_ids);
plot_recon_batch(fullfile(DIAG_DIR, 'avg_recon'),   'Avg',   avg_k,   per_sample_rel_val, val_idx, U3, U3_val_recon, s_grid, sample_ids);

fprintf('  Saved %d worst, %d best, %d avg plots\n', numel(worst_k), numel(best_k), numel(avg_k));

diag.q5_worst_ids = sample_ids(val_idx(worst_k))';
diag.q5_worst_rel = per_sample_rel_val(worst_k)';
diag.q5_best_ids  = sample_ids(val_idx(best_k))';
diag.q5_best_rel  = per_sample_rel_val(best_k)';
diag.q5_avg_ids   = sample_ids(val_idx(avg_k))';
diag.q5_avg_rel   = per_sample_rel_val(avg_k)';

%% =========================================================================
%% Q6 — Learnability baseline (ridge regression from cheap inputs)
%% =========================================================================
fprintf('\n--- Q6: Learnability baseline (ridge on coarse inputs) ---\n');

ang_train = t.ang_deg(train_idx);
ang_val   = t.ang_deg(val_idx);
tr_train  = t.tr_ratio(train_idx);
tr_val    = t.tr_ratio(val_idx);

X_train = [tr_train, ang_train, sin(ang_train*pi/180), cos(ang_train*pi/180)];
X_val   = [tr_val,   ang_val,   sin(ang_val*pi/180),   cos(ang_val*pi/180)];

n_check = min(4, K);
r2_val = zeros(n_check, 1);
for k = 1:n_check
    mdl = fitrlinear(X_train, Z_train(:,k), 'Learner','leastsquares', 'Regularization','ridge', 'Lambda', 1e-3);
    y_pred = predict(mdl, X_val);
    ss_res = sum((Z_val(:,k) - y_pred).^2);
    ss_tot = sum((Z_val(:,k) - mean(Z_train(:,k))).^2);
    r2_val(k) = 1 - ss_res/ss_tot;
    fprintf('  PC%d  ridge R² (val) = %.4f\n', k, r2_val(k));
end

q6_pass = r2_val(1) > 0;
fprintf('  Q6: %s\n', verdict(q6_pass));

for k = 1:n_check
    diag.(sprintf('q6_ridge_r2_pc%d_val', k)) = r2_val(k);
end
diag.q6_pass = q6_pass;

%% =========================================================================
%% Summary + manifest
%% =========================================================================
fprintf('\n========= SUMMARY =========\n');
fprintf('  Q1 reconstruction fidelity : %s\n', verdict(q1_pass));
fprintf('  Q2 variance compactness    : %s\n', verdict(q2_pass));
fprintf('  Q3 basis stability         : %s\n', verdict(q3_pass));
fprintf('  Q6 learnability baseline   : %s\n', verdict(q6_pass));
fprintf('===========================\n\n');

all_pass = q1_pass && q2_pass && q3_pass && q6_pass;
if all_pass
    fprintf('All checks PASS — PCA target is usable.\n');
else
    fprintf('One or more checks WARN — review diagnostics before training.\n');
end

try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

diag.n_components        = K;
diag.all_pass            = all_pass;
diag.git_sha             = git_sha;
diag.created_at          = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

diagFile = fullfile(DIAG_DIR, 'pca_diagnostics.json');
fid = fopen(diagFile, 'w');
fprintf(fid, '%s\n', jsonencode(diag));
fclose(fid);
fprintf('Diagnostics written: %s\n', diagFile);

%% ---- Local helpers --------------------------------------------------------
function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'WARN'; end
end

function plot_recon_batch(out_dir, label, indices_in_val, per_sample_rel_val, val_idx, U3, U3_val_recon, s_grid, sample_ids)
if ~isfolder(out_dir), mkdir(out_dir); end
for wi = 1:numel(indices_in_val)
    idx_in_val = indices_in_val(wi);
    global_idx = val_idx(idx_in_val);
    sid = sample_ids{global_idx};
    fig = figure('Visible','off');
    fig.Position = [0 0 600 420];
    plot(s_grid, U3(global_idx,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'True');
    hold on;
    plot(s_grid, U3_val_recon(idx_in_val,:), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Recon');
    title(sprintf('%s rank %d  —  Rel err = %.2f%%\n%s', label, wi, ...
        100*per_sample_rel_val(idx_in_val), strrep(sid,'_','\_')), ...
        'FontSize', 9, 'Interpreter', 'tex');
    xlabel('s (normalized arc length)'); ylabel('u3 (m)');
    legend('Location','best');
    grid on;
    fname = fullfile(out_dir, sprintf('%s_%02d_%s.png', lower(label), wi, sid));
    exportgraphics(fig, fname, 'Resolution', 150);
    close(fig);
end
end
