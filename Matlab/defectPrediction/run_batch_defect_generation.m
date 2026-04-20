% RUN_BATCH_DEFECT_GENERATION Generate defect-only PSSCone cases.
%
% How to run:
%   run('Matlab/defectPrediction/run_batch_defect_generation.m')
%
% Results are written to:
%   Matlab/defectPrediction/results/lam{ang}/def{frac}/

cfg = defect_prediction_config();
addpath(cfg.matlabRoot);

if ~isfolder(cfg.resultsRoot)
    mkdir(cfg.resultsRoot);
end

if ~isfolder(cfg.stagingRoot)
    mkdir(cfg.stagingRoot);
end

fprintf('Generating %d cases (angles x fractions) under %s\n', ...
    numel(cfg.cases), cfg.resultsRoot);

for i = 1:numel(cfg.cases)
    c = cfg.cases(i);

    if isfolder(c.case_dir)
        fprintf('[SKIP] %s already exists: %s\n', c.case_name, c.case_dir);
        continue;
    end

    % Ensure lam{ang} parent dir exists
    angDir = fileparts(c.case_dir);
    if ~isfolder(angDir)
        mkdir(angDir);
    end

    stageRoot = tempname(cfg.stagingRoot);
    mkdir(stageRoot);

    params = cfg.baseParams;
    params.lamellarAngleDeg    = c.lamellar_angle_deg;
    params.remove_top_spin_frac = c.defect_fraction;
    params.resultsRoot          = stageRoot;

    fprintf('[RUN ] %s | angle=%d deg | remove_top_spin_frac=%.1f\n', ...
        c.case_name, c.lamellar_angle_deg, c.defect_fraction);
    PSSCone(params);

    generatedRunDir = find_single_generated_run(stageRoot);
    originalRunFolder = get_last_path_part(generatedRunDir);
    move_ok = movefile(generatedRunDir, c.case_dir);
    if ~move_ok
        error('run_batch_defect_generation:MoveFailed', ...
            'Failed to move %s to %s.', generatedRunDir, c.case_dir);
    end

    write_case_metadata(c.case_dir, c.case_name, c.defect_fraction, ...
        c.lamellar_angle_deg, originalRunFolder);
    fprintf('[DONE] %s -> %s\n', c.case_name, c.case_dir);
end

function runDir = find_single_generated_run(stageRoot)
entries = dir(stageRoot);
entries = entries([entries.isdir]);
entries = entries(~ismember({entries.name}, {'.', '..'}));

if numel(entries) ~= 1
    error('run_batch_defect_generation:UnexpectedOutput', ...
        'Expected exactly one generated run folder under %s, found %d.', ...
        stageRoot, numel(entries));
end

runDir = fullfile(stageRoot, entries(1).name);
end

function write_case_metadata(caseDir, caseName, defectFrac, angleDeg, originalRunFolder)
meta = struct();
meta.case_name            = caseName;
meta.remove_top_spin_frac = defectFrac;
meta.lamellar_angle_deg   = angleDeg;
meta.original_run_folder  = originalRunFolder;
meta.generator            = 'PSSCone';
meta.generated_at         = datestr(now, 31);

metaPath = fullfile(caseDir, 'defect_case.json');
fid = fopen(metaPath, 'w');
if fid < 0
    error('run_batch_defect_generation:MetaWriteFailed', ...
        'Unable to write case metadata to %s.', metaPath);
end

cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(meta));
end

function name = get_last_path_part(pathStr)
[~, name] = fileparts(pathStr);
end
