% step5_build_u3_target_dataset.m
%
% Builds a per-sample u3(s) midpoint displacement profile for every sample
% in data/dataset/samples/ and aggregates QC-passed profiles into
% data/dataset/targets/u3_targets.mat.
%
% Inputs  (per sample):  midpoint_results_shell.csv
% Outputs (per sample):  target_u3_profile.mat
% Outputs (aggregated):  data/dataset/targets/u3_targets.mat
%                        data/dataset/targets/u3_targets_manifest.json
% Log:                   logs/step5_u3_target_log.txt
%
% Target design:
%   Arc length is computed on the DEFORMED curve (Y+U2, Z+U3).
%   Undeformed Y is used only to obtain a stable node ordering.
%   u3(s) and u2(s) are resampled on a fixed s_grid in [0,1].
%   PCA is NOT fit here; U3 matrix rows are QC-passed-only so a later
%   step6 can slice a stratified train split and call pca(U3(train_idx,:)).

%% ---- Config ---------------------------------------------------------------
REPO_ROOT    = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
DATASET_DIR  = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'samples');
TARGETS_DIR  = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'targets');
LOG_FILE     = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'logs', 'step5_u3_target_log.txt');

N_GRID          = 128;       % resampling resolution
MIN_RAW_PTS     = 64;        % minimum raw midplane points
U3_MAX_ABS_M    = 0.5;       % sanity cap: |u3| > 0.5 m triggers fail:u3_out_of_range
ZERO_U_TOL      = 1e-14;     % all-zero displacement + stress means failed/no-result extraction
ZERO_STRESS_TOL = 1e-6;
SCHEMA_VERSION  = 1;
OVERWRITE       = false;     % set true to reprocess already-completed samples

%% ---- Init -----------------------------------------------------------------
ensure_dir(TARGETS_DIR);
ensure_dir(fileparts(LOG_FILE));

logFid = fopen(LOG_FILE, 'a');
log_entry(logFid, 'INFO', 'step5 started', '', '', '');

created_at = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

%% ---- Discover sample folders ----------------------------------------------
items = dir(DATASET_DIR);
sampleDirs = {};
for k = 1:numel(items)
    if items(k).isdir && items(k).name(1) ~= '.'
        sampleDirs{end+1,1} = fullfile(DATASET_DIR, items(k).name); %#ok<AGROW>
    end
end
n_total = numel(sampleDirs);
fprintf('Found %d sample folders in %s\n', n_total, DATASET_DIR);
if n_total == 0
    error('step5:no_samples', 'No sample folders under %s', DATASET_DIR);
end

%% ---- Per-sample processing (serialized; parfor-ready if needed) -----------
fprintf('[1/3] Processing %d samples...\n', n_total);
t_loop_start = tic;
n_skipped = 0; n_processed = 0; n_failed = 0;
PRINT_EVERY = 100;

default_result = struct( ...
    'sample_id',  '', ...
    'tr_ratio',   NaN, ...
    'ang_deg',    NaN, ...
    'run_suffix', '', ...
    'qc_status',  '', ...
    'u3',         [], ...
    'u2',         [], ...
    'ok',         false);
results = repmat(default_result, n_total, 1);

for i = 1:n_total
    sDir = sampleDirs{i};
    [~, folderName] = fileparts(sDir);

    % Folder name is the canonical sample id. This avoids thousands of
    % optional metadata reads on synced/cloud folders.
    sample_id = folderName;
    results(i).sample_id = sample_id;

    % --- parse stratification keys from sample_id ---------------------------
    [tr_ratio, ang_deg, run_suffix] = parse_sample_id(sample_id);
    results(i).tr_ratio   = tr_ratio;
    results(i).ang_deg    = ang_deg;
    results(i).run_suffix = run_suffix;

    % --- skip if already done ------------------------------------------------
    outMat = fullfile(sDir, 'target_u3_profile.mat');
    if ~OVERWRITE && isfile(outMat)
        % reload qc_status to include in aggregate
        try
            d = load(outMat, 'qc_status', 'u3', 'u2');
            results(i).qc_status = d.qc_status;
            results(i).ok        = strcmp(d.qc_status, 'ok');
            if results(i).ok
                results(i).u3 = d.u3;
                results(i).u2 = d.u2;
            end
        catch
            results(i).qc_status = 'fail:load_existing_error';
        end
        n_skipped = n_skipped + 1;
    if mod(i, PRINT_EVERY) == 0
        fprintf('  [1/3] %d/%d  (skipped=%d processed=%d failed=%d)  %.0fs elapsed\n', ...
            i, n_total, n_skipped, n_processed, n_failed, toc(t_loop_start));
        drawnow('limitrate');
    end
        continue
    end

    % --- load CSV ------------------------------------------------------------
    csvPath = fullfile(sDir, 'midpoint_results_shell.csv');
    if ~isfile(csvPath)
        qc = 'fail:missing_csv';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    try
        M = readmatrix(csvPath, 'NumHeaderLines', 1);
    catch err
        qc = 'fail:csv_parse_error';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, err.message);
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- extract columns -----------------------------------------------------
    try
        if size(M, 2) < 7
            error('expected at least 7 numeric columns, got %d', size(M, 2));
        end
        Y  = M(:, 3);
        Z  = M(:, 4);
        U2 = M(:, 6);
        U3 = M(:, 7);
        if size(M, 2) >= 14
            S_Mises = M(:, 14);
        else
            S_Mises = [];
        end
    catch
        qc = 'fail:missing_columns';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, 'Expected Y,Z,U2,U3 columns');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- QC: NaN / Inf -------------------------------------------------------
    if any(~isfinite(Y)) || any(~isfinite(Z)) || any(~isfinite(U2)) || any(~isfinite(U3))
        qc = 'fail:nan_inf';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- QC: failed/no-result Abaqus extraction -------------------------------
    max_u = max([abs(U2); abs(U3)]);
    max_stress = 0;
    if ~isempty(S_Mises)
        max_stress = max(abs(S_Mises));
    end
    if max_u <= ZERO_U_TOL && max_stress <= ZERO_STRESS_TOL
        qc = 'fail:zero_displacement_and_stress';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- QC: enough points ---------------------------------------------------
    n_pts = numel(Y);
    if n_pts < MIN_RAW_PTS
        qc = sprintf('fail:too_few_points(%d)', n_pts);
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- QC: u3 magnitude sanity check ---------------------------------------
    if max(abs(U3)) > U3_MAX_ABS_M
        qc = sprintf('fail:u3_out_of_range(%.4f m)', max(abs(U3)));
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        n_failed = n_failed + 1;
        save_per_sample_fail(outMat, sample_id, qc, tr_ratio, ang_deg, run_suffix, SCHEMA_VERSION, created_at);
        continue
    end

    % --- Order by undeformed Y (stable node ordering) -----------------------
    [~, ord] = sort(Y);
    y_raw  = Y(ord);
    z_raw  = Z(ord);
    u2_raw = U2(ord);
    u3_raw = U3(ord);

    % --- Deformed arc length on (Y+U2, Z+U3) --------------------------------
    yd = y_raw + u2_raw;
    zd = z_raw + u3_raw;
    ds = sqrt(diff(yd).^2 + diff(zd).^2);
    arc_cumul = [0; cumsum(ds)];
    arc_length_m = arc_cumul(end);

    % --- QC: monotone cumulative arc length ----------------------------------
    if any(diff(arc_cumul) <= 0)
        qc = 'fail:non_monotone_arc_length';
        log_entry(logFid, 'FAIL', sample_id, qc, csvPath, '');
        results(i).qc_status = qc;
        s_raw = arc_cumul / arc_length_m;
        n_failed = n_failed + 1;
        save_per_sample_with_raw(outMat, sample_id, qc, ...
            y_raw, z_raw, u2_raw, u3_raw, s_raw, arc_length_m, n_pts, ...
            [], [], [], ...
            tr_ratio, ang_deg, run_suffix, N_GRID, SCHEMA_VERSION, created_at);
        continue
    end

    % --- Normalize arc length ------------------------------------------------
    s_raw = arc_cumul / arc_length_m;

    % --- Resample on fixed grid ----------------------------------------------
    s_grid = linspace(0, 1, N_GRID);
    u3_res = interp1(s_raw, u3_raw, s_grid, 'pchip');
    u2_res = interp1(s_raw, u2_raw, s_grid, 'pchip');

    % --- Write per-sample output --------------------------------------------
    qc = 'ok';
    save_per_sample_with_raw(outMat, sample_id, qc, ...
        y_raw, z_raw, u2_raw, u3_raw, s_raw, arc_length_m, n_pts, ...
        s_grid, u3_res, u2_res, ...
        tr_ratio, ang_deg, run_suffix, N_GRID, SCHEMA_VERSION, created_at);

    results(i).qc_status = qc;
    results(i).ok        = true;
    results(i).u3        = u3_res;
    results(i).u2        = u2_res;
    n_processed = n_processed + 1;
        if mod(i, PRINT_EVERY) == 0
            fprintf('  [1/3] %d/%d  (skipped=%d processed=%d failed=%d)  %.0fs elapsed\n', ...
                i, n_total, n_skipped, n_processed, n_failed, toc(t_loop_start));
            drawnow('limitrate');
        end
end
fprintf('[1/3] Loop done in %.1fs — skipped=%d  ok=%d  failed=%d\n', ...
    toc(t_loop_start), n_skipped, n_processed, n_failed);

%% ---- Aggregate QC-passed samples only -------------------------------------
fprintf('[2/3] Aggregating QC-passed samples into u3_targets.mat...\n');
schema_version = SCHEMA_VERSION;
s_grid = linspace(0, 1, N_GRID);

ok_idx = find([results.ok]);
M = numel(ok_idx);

sample_ids = {results(ok_idx).sample_id}';
tr_ratio   = [results(ok_idx).tr_ratio]';
ang_deg    = [results(ok_idx).ang_deg]';
run_suffix = {results(ok_idx).run_suffix}';

U3_mat = zeros(M, N_GRID);
U2_mat = zeros(M, N_GRID);
for j = 1:M
    U3_mat(j,:) = results(ok_idx(j)).u3;
    U2_mat(j,:) = results(ok_idx(j)).u2;
end

targetsFile = fullfile(TARGETS_DIR, 'u3_targets.mat');
save(targetsFile, 'sample_ids', 'U3_mat', 'U2_mat', 's_grid', ...
    'tr_ratio', 'ang_deg', 'run_suffix', ...
    'schema_version', 'created_at', '-v7.3'); %#ok<USENS>
fprintf('[2/3] Saved aggregated target: %s  (%d/%d QC-passed)\n', targetsFile, M, n_total);

%% ---- Manifest -------------------------------------------------------------
all_statuses = {results.qc_status};
fail_reasons = {};
for i = 1:numel(all_statuses)
    s = all_statuses{i};
    if ~strcmp(s, 'ok')
        fail_reasons{end+1} = s; %#ok<AGROW>
    end
end

unique_reasons = unique(fail_reasons);
n_failed_by_reason = struct();
for k = 1:numel(unique_reasons)
    key = matlab_safe_fieldname(unique_reasons{k});
    n_failed_by_reason.(key) = sum(strcmp(fail_reasons, unique_reasons{k}));
end

try
    [~, git_sha_raw] = system('git rev-parse HEAD');
    git_sha = strtrim(git_sha_raw);
catch
    git_sha = 'unknown';
end

manifest = struct( ...
    'n_samples_total',   n_total, ...
    'n_ok',              M, ...
    'n_failed',          n_total - M, ...
    'n_failed_by_reason', n_failed_by_reason, ...
    'n_grid',            N_GRID, ...
    'schema_version',    SCHEMA_VERSION, ...
    'min_raw_pts',       MIN_RAW_PTS, ...
    'u3_max_abs_m',      U3_MAX_ABS_M, ...
    'zero_u_tol',        ZERO_U_TOL, ...
    'zero_stress_tol',   ZERO_STRESS_TOL, ...
    'arc_length_basis',  'deformed_curve_YpU2_ZpU3', ...
    'git_sha',           git_sha, ...
    'created_at',        created_at);

fprintf('[3/3] Writing manifest...\n');
manifestFile = fullfile(TARGETS_DIR, 'u3_targets_manifest.json');
fid = fopen(manifestFile, 'w');
fprintf(fid, '%s\n', jsonencode(manifest));
fclose(fid);
fprintf('[3/3] Manifest written: %s\n', manifestFile);

log_entry(logFid, 'INFO', 'step5 complete', ...
    sprintf('ok=%d fail=%d total=%d', M, n_total-M, n_total), '', '');
fclose(logFid);
fprintf('Done. %d/%d samples passed QC.\n', M, n_total);

%% =========================================================================
%% Local helpers
%% =========================================================================

function save_per_sample_fail(outMat, sample_id, qc_status, ...
        tr_ratio, ang_deg, run_suffix, schema_version, created_at)
save(outMat, 'sample_id', 'qc_status', ...
    'tr_ratio', 'ang_deg', 'run_suffix', ...
    'schema_version', 'created_at');
end

function save_per_sample_with_raw(outMat, sample_id, qc_status, ...
        y_raw, z_raw, u2_raw, u3_raw, s_raw, arc_length_m, n_points_raw, ...
        s_grid, u3, u2, ...
        tr_ratio, ang_deg, run_suffix, N_GRID, schema_version, created_at) %#ok<INUSL>
save(outMat, ...
    'sample_id', 'qc_status', ...
    's_grid', 'u3', 'u2', ...
    'n_points_raw', 'arc_length_m', ...
    'tr_ratio', 'ang_deg', 'run_suffix', ...
    'schema_version', 'created_at');
end

function [tr_ratio, ang_deg, run_suffix] = parse_sample_id(sample_id)
% Parse stratification keys from a run_name like:
%   sheetCone_tr100_ang045_lamellar_N128_1x1
%   sheetCone_tr050_ang090_lamellar_N128_1x1_run03
tr_ratio   = NaN;
ang_deg    = NaN;
run_suffix = 'run01';

tk = regexp(sample_id, 'tr(\d+)', 'tokens', 'once');
if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end

tk = regexp(sample_id, 'ang(\d+)', 'tokens', 'once');
if ~isempty(tk), ang_deg = str2double(tk{1}); end

tk = regexp(sample_id, '_run(\d+)$', 'tokens', 'once');
if ~isempty(tk), run_suffix = ['run' tk{1}]; end
end

function key = matlab_safe_fieldname(str)
key = regexprep(str, '[^a-zA-Z0-9_]', '_');
if ~isempty(key) && ~isnan(str2double(key(1)))
    key = ['x' key];
end
end

function ensure_dir(d)
if ~isfolder(d), mkdir(d); end
end

function log_entry(fid, level, sample_id, reason, artifact_path, detail)
fprintf(fid, '[%s] %s | sample=%s | reason=%s | path=%s | detail=%s\n', ...
    datestr(now,'yyyy-mm-ddTHH:MM:SS'), level, sample_id, reason, artifact_path, detail);
end
