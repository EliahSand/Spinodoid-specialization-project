% RUN_BATCH_DEFECT_GENERATION Generate defect-only PSSCone cases.
%
% How to run:
%   run('Matlab/defectPrediction/run_batch_defect_generation.m')
%

cfg = defect_prediction_config();
addpath(cfg.matlabRoot);

if ~isfolder(cfg.resultsRoot)
    mkdir(cfg.resultsRoot);
end

if ~isfolder(cfg.stagingRoot)
    mkdir(cfg.stagingRoot);
end

fprintf('Generating %d defect cases under %s\n', numel(cfg.caseNames), cfg.resultsRoot);

for i = 1:numel(cfg.defectFractions)
    caseName = cfg.caseNames{i};
    caseDir = cfg.caseDirs{i};
    defectFrac = cfg.defectFractions(i);

    if isfolder(caseDir)
        fprintf('[SKIP] %s already exists: %s\n', caseName, caseDir);
        continue;
    end

    stageRoot = tempname(cfg.stagingRoot);
    mkdir(stageRoot);

    params = cfg.baseParams;
    params.remove_top_spin_frac = defectFrac;
    params.resultsRoot = stageRoot;

    fprintf('[RUN ] %s | remove_top_spin_frac=%.1f\n', caseName, defectFrac);
    PSSCone(params);

    generatedRunDir = find_single_generated_run(stageRoot);
    originalRunFolder = get_last_path_part(generatedRunDir);
    move_ok = movefile(generatedRunDir, caseDir);
    if ~move_ok
        error('run_batch_defect_generation:MoveFailed', ...
            'Failed to move %s to %s.', generatedRunDir, caseDir);
    end

    write_case_metadata(caseDir, caseName, defectFrac, originalRunFolder);
    fprintf('[DONE] %s -> %s\n', caseName, caseDir);
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

function write_case_metadata(caseDir, caseName, defectFrac, originalRunFolder)
meta = struct();
meta.case_name = caseName;
meta.remove_top_spin_frac = defectFrac;
meta.original_run_folder = originalRunFolder;
meta.generator = 'PSSCone';
meta.generated_at = datestr(now, 31);

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
