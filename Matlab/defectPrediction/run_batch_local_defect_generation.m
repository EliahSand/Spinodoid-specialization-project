% RUN_BATCH_LOCAL_DEFECT_GENERATION  Generate the 35-case local-defect sweep.
%
%   Runs PSSCone for each (angle, defect_type, position) case, then stamps
%   cracks or holes into the top spin layer in a post-processing step.
%   No GRF removal is applied (remove_top_spin_frac = 0 always).
%
% How to run:
%   run('Matlab/defectPrediction/run_batch_local_defect_generation.m')
%
% Results under:
%   Matlab/defectPrediction/results/lam{ang}/{defect_type}_{position}/

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);  % apply_local_defects, export_sheet_stl

cfg = defect_prediction_config();
addpath(cfg.matlabRoot);

if ~isfolder(cfg.resultsRoot), mkdir(cfg.resultsRoot); end
if ~isfolder(cfg.stagingRoot), mkdir(cfg.stagingRoot); end

cases = cfg.local_cases;
fprintf('Generating %d local-defect cases under %s\n', numel(cases), cfg.resultsRoot);

for i = 1:numel(cases)
    c = cases(i);

    if isfolder(c.case_dir)
        fprintf('[SKIP] %s\n', c.case_name);
        continue;
    end

    angDir = fileparts(c.case_dir);
    if ~isfolder(angDir), mkdir(angDir); end

    stageRoot = tempname(cfg.stagingRoot);
    mkdir(stageRoot);

    % Run PSSCone with GRF disabled.
    params = cfg.baseParams;
    params.lamellarAngleDeg    = c.lamellar_angle_deg;
    params.remove_top_spin_frac = 0.0;
    params.resultsRoot          = stageRoot;

    fprintf('[RUN ] %s | angle=%d | type=%s | pos=%s\n', ...
        c.case_name, c.lamellar_angle_deg, c.defect_type, c.position_name);
    PSSCone(params);

    generatedRunDir = local_find_single_run(stageRoot);

    % --- Post-process: stamp local defect into the spin layer. ---
    matPath = fullfile(generatedRunDir, 'sheet.mat');
    S = load(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness');
    sheetMask       = S.sheetMask;
    voxelSizeXY     = S.voxelSizeXY;
    zVoxelThickness = S.zVoxelThickness;

    L    = params.L;
    N    = params.N;
    Svox = L / N;
    tbV  = max(3, round(params.t_base / Svox));
    tsV  = max(1, round(params.t_spin / Svox));

    spinLayer = sheetMask(:,:,tbV+1:end);

    if ~isempty(c.local_defects)
        spinLayer = apply_local_defects(spinLayer, Svox, L, c.local_defects);
    end

    sheetMask(:,:,tbV+1:end) = spinLayer;

    % Verify periodicity.
    xMismatch = nnz(sheetMask(1,:,:) ~= sheetMask(end,:,:));
    yMismatch = nnz(sheetMask(:,1,:) ~= sheetMask(:,end,:));
    if xMismatch || yMismatch
        error('run_batch_local_defect_generation:periodicity', ...
            '%s: periodicity broken after defect stamp (x=%d, y=%d).', ...
            c.case_name, xMismatch, yMismatch);
    end

    save(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness', '-v7.3');

    % Export STL.
    Lz = params.t_base + params.t_spin;
    export_sheet_stl(generatedRunDir, sheetMask, Svox, Lz);

    % Move to final location.
    ok = movefile(generatedRunDir, c.case_dir);
    if ~ok
        error('run_batch_local_defect_generation:MoveFailed', ...
            'Failed to move %s to %s.', generatedRunDir, c.case_dir);
    end

    local_write_metadata(c);
    fprintf('[DONE] %s -> %s\n', c.case_name, c.case_dir);
end

% -----------------------------------------------------------------------
function runDir = local_find_single_run(stageRoot)
entries = dir(stageRoot);
entries = entries([entries.isdir]);
entries = entries(~ismember({entries.name}, {'.', '..'}));
if numel(entries) ~= 1
    error('run_batch_local_defect_generation:UnexpectedOutput', ...
        'Expected 1 run folder under %s, found %d.', stageRoot, numel(entries));
end
runDir = fullfile(stageRoot, entries(1).name);
end

function local_write_metadata(c)
meta = struct();
meta.case_name          = c.case_name;
meta.lamellar_angle_deg = c.lamellar_angle_deg;
meta.defect_type        = c.defect_type;
meta.position_name      = c.position_name;
meta.position_xy        = c.position_xy;
meta.remove_top_spin_frac = 0.0;
meta.generator          = 'PSSCone+local';
meta.generated_at       = datestr(now, 31);

if ~isempty(c.severity) && ~isempty(fieldnames(c.severity))
    meta.severity = c.severity;
end
if ~isempty(c.local_defects)
    meta.local_defects = c.local_defects;
end

metaPath = fullfile(c.case_dir, 'defect_case.json');
fid = fopen(metaPath, 'w');
if fid < 0
    error('run_batch_local_defect_generation:MetaFailed', ...
        'Cannot write metadata to %s.', metaPath);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(meta));
end
