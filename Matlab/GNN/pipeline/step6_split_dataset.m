% step6_split_dataset.m
%
% Builds a reproducible stratified train/val/test split over the QC-passed
% samples in u3_targets.mat and saves indices to split_indices.mat.
%
% Stratification: each (tr_ratio, ang_deg) bin is split independently so
% that all thickness ratios and angles are represented in every partition.
% This is necessary because samples are stored in increasing lamellar order
% — a contiguous split would produce a severely biased dataset.
%
% Inputs:   data/dataset/targets/u3_targets.mat
% Outputs:  data/dataset/targets/split_indices.mat
%           data/dataset/targets/split_manifest.json

%% ---- Config ---------------------------------------------------------------
REPO_ROOT   = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'targets');

TRAIN_FRAC = 0.70;
VAL_FRAC   = 0.15;
% TEST_FRAC  = 1 - TRAIN_FRAC - VAL_FRAC = 0.15

RNG_SEED = 42;

%% ---- Load -----------------------------------------------------------------
targetsFile = fullfile(TARGETS_DIR, 'u3_targets.mat');
if ~isfile(targetsFile)
    error('step6:missing_targets', 'Run step5 first: %s not found.', targetsFile);
end

fprintf('Loading %s...\n', targetsFile);
d = load(targetsFile, 'sample_ids', 'tr_ratio', 'ang_deg');
M = numel(d.sample_ids);
fprintf('Loaded %d QC-passed samples.\n', M);

%% ---- Stratified split -----------------------------------------------------
rng(RNG_SEED);

train_idx = [];
val_idx   = [];
test_idx  = [];

tr_vals  = unique(d.tr_ratio);
ang_vals = unique(d.ang_deg);

for ti = 1:numel(tr_vals)
    for ai = 1:numel(ang_vals)
        bin_mask = d.tr_ratio == tr_vals(ti) & d.ang_deg == ang_vals(ai);
        bin_idx  = find(bin_mask);
        if isempty(bin_idx), continue; end

        % shuffle within bin
        bin_idx = bin_idx(randperm(numel(bin_idx)));
        n = numel(bin_idx);

        n_train = max(1, round(n * TRAIN_FRAC));
        n_val   = max(1, round(n * VAL_FRAC));
        n_test  = n - n_train - n_val;

        if n_test < 1
            % bin too small to fill all three splits: give test priority
            n_train = max(1, n - 2);
            n_val   = 1;
            n_test  = max(1, n - n_train - n_val);
            n_train = n - n_val - n_test;
        end

        train_idx = [train_idx; bin_idx(1:n_train)];                     %#ok<AGROW>
        val_idx   = [val_idx;   bin_idx(n_train+1:n_train+n_val)];       %#ok<AGROW>
        test_idx  = [test_idx;  bin_idx(n_train+n_val+1:end)];           %#ok<AGROW>
    end
end

% sort so indices are in a deterministic order
train_idx = sort(train_idx);
val_idx   = sort(val_idx);
test_idx  = sort(test_idx);

n_train = numel(train_idx);
n_val   = numel(val_idx);
n_test  = numel(test_idx);

fprintf('Split: train=%d (%.1f%%)  val=%d (%.1f%%)  test=%d (%.1f%%)\n', ...
    n_train, 100*n_train/M, n_val, 100*n_val/M, n_test, 100*n_test/M);

assert(n_train + n_val + n_test == M, 'Split indices do not cover all samples.');
assert(numel(unique([train_idx; val_idx; test_idx])) == M, 'Duplicate indices in split.');

%% ---- Save -----------------------------------------------------------------
created_at    = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
schema_version = 1;
rng_seed       = RNG_SEED;
train_frac     = TRAIN_FRAC;
val_frac       = VAL_FRAC;
test_frac      = 1 - TRAIN_FRAC - VAL_FRAC;

splitFile = fullfile(TARGETS_DIR, 'split_indices.mat');
save(splitFile, 'train_idx', 'val_idx', 'test_idx', ...
    'rng_seed', 'train_frac', 'val_frac', 'test_frac', ...
    'schema_version', 'created_at');
fprintf('Saved split indices: %s\n', splitFile);

%% ---- Manifest -------------------------------------------------------------
try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

manifest = struct( ...
    'n_total',        M, ...
    'n_train',        n_train, ...
    'n_val',          n_val, ...
    'n_test',         n_test, ...
    'train_frac_actual', n_train/M, ...
    'val_frac_actual',   n_val/M, ...
    'test_frac_actual',  n_test/M, ...
    'rng_seed',       RNG_SEED, ...
    'stratify_by',    'tr_ratio x ang_deg', ...
    'schema_version', schema_version, ...
    'git_sha',        git_sha, ...
    'created_at',     created_at);

manifestFile = fullfile(TARGETS_DIR, 'split_manifest.json');
fid = fopen(manifestFile, 'w');
fprintf(fid, '%s\n', jsonencode(manifest));
fclose(fid);
fprintf('Manifest: %s\n', manifestFile);
fprintf('Done.\n');
