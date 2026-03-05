function summary = run_all_analysis(opts)
%RUN_ALL_ANALYSIS Run overlay + curvature analysis for all detected runs.
%
% Finds run folders under results/sheets that contain both FEA and
% FEA_shell CSV outputs, then calls compare_single_run for each.

arguments
    opts.ResultsSheetsRoot (1, :) char = ''
    opts.OutputRoot (1, :) char = ''
    opts.DoPlots (1, 1) logical = true
    opts.WriteFiles (1, 1) logical = true
    opts.MakeOverlaySummary (1, 1) logical = true
end

tStart = tic;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');

if isempty(opts.ResultsSheetsRoot)
    resultsSheetsRoot = fullfile(scriptDir, '..', '..', 'results', 'sheets');
else
    resultsSheetsRoot = resolve_path(opts.ResultsSheetsRoot, repoRoot);
end
resultsSheetsRoot = char(resultsSheetsRoot);

if isempty(opts.OutputRoot)
    outputRoot = fullfile(scriptDir, '..', '..', 'results', 'analysis', 'single_runs');
else
    outputRoot = resolve_path(opts.OutputRoot, repoRoot);
end
outputRoot = char(outputRoot);

ensure_dir(outputRoot);
runFolders = discover_run_folders(resultsSheetsRoot);

if isempty(runFolders)
    warning('run_all_analysis:NoRunsFound', ...
        'No run folders with both FEA and FEA_shell CSV files were found under %s', ...
        resultsSheetsRoot);
    summary = table();
    return;
end

fprintf('Found %d run folder(s). Starting analysis...\n', numel(runFolders));

overlayFields = {'U1','U2','U3','S_11','S_22','S_33','S_12','S_13','S_23','S_Mises'};
overlayData = init_overlay_data(overlayFields);

rows = {};
for i = 1:numel(runFolders)
    runFolder = runFolders{i};
    [~, runName] = fileparts(runFolder);
    runOutputDir = fullfile(outputRoot, runName);
    fprintf('\n[%d/%d] %s\n', i, numel(runFolders), runName);

    status = 'ok';
    kTheory = NaN;
    kSolid = NaN;
    kShell = NaN;
    relSolid = NaN;
    relShell = NaN;
    trLabel = '';
    thetaDeg = NaN;

    try
        compare_single_run(runFolder, ...
            'OutputDir', runOutputDir, ...
            'DoPlots', opts.DoPlots, ...
            'WriteFiles', opts.WriteFiles);

        [solid, shell] = load_run_pair_local(runFolder);
        [trLabelChar, thetaDeg] = parse_run_tags_local(runName, runFolder);
        trLabel = trLabelChar;

        [curv, ~] = check_curvature_single_run(runFolder, solid, shell);
        kTheory = curv.theory.kappa;
        kSolid = curv.solid.kappaCircle;
        kShell = curv.shell.kappaCircle;
        relSolid = curv.solid.RelErrTheory;
        relShell = curv.shell.RelErrTheory;

        overlayData = append_overlay_data(overlayData, overlayFields, trLabelChar, thetaDeg, solid, shell);
    catch ME
        status = 'failed';
        warning('run_all_analysis:RunFailed', ...
            'Run failed for %s:\n%s', runFolder, ...
            ME.getReport('extended', 'hyperlinks', 'off'));
    end

    rows(end+1, :) = {runName, runFolder, trLabel, thetaDeg, status, ...
        kTheory, kSolid, kShell, relSolid, relShell}; %#ok<AGROW>
end

summary = cell2table(rows, 'VariableNames', ...
    {'RunName','RunFolder','RatioTag','ThetaDeg','Status', ...
     'KappaTheory_1_per_m','KappaSolid_1_per_m','KappaShell_1_per_m', ...
     'RelErrSolidToTheory','RelErrShellToTheory'});

writetable(summary, fullfile(outputRoot, 'run_all_summary.csv'));

if opts.MakeOverlaySummary && opts.DoPlots
    overlayOutDir = fullfile(outputRoot, '_overlay_summary');
    plot_overlay_series(overlayData, overlayOutDir);
    fprintf('Saved overlay summary plots to: %s\n', overlayOutDir);
end

fprintf('\nFinished run_all_analysis in %.2f s.\n', toc(tStart));
end

function pOut = resolve_path(pIn, repoRoot)
pOut = pIn;
if isempty(pIn)
    return;
end
if is_absolute_path(pIn)
    return;
end
pOut = fullfile(repoRoot, pIn);
end

function tf = is_absolute_path(p)
if isempty(p)
    tf = false;
    return;
end
if ispc
    tf = ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')) || startsWith(p, '\\');
else
    tf = startsWith(p, '/');
end
end

function runFolders = discover_run_folders(resultsSheetsRoot)
runFolders = {};
if ~isfolder(resultsSheetsRoot)
    return;
end

tmp = {};
allDirs = collect_subdirs(resultsSheetsRoot);
for i = 1:numel(allDirs)
    runFolder = allDirs{i};
    shellDir = fullfile(runFolder, 'FEA_shell');
    if ~isfolder(shellDir)
        continue;
    end
    solidCsv = dir(fullfile(runFolder, 'FEA', '*.csv'));
    shellCsv = dir(fullfile(shellDir, '*.csv'));
    if isempty(solidCsv) || isempty(shellCsv)
        continue;
    end
    tmp{end+1,1} = runFolder; %#ok<AGROW>
end

if isempty(tmp)
    runFolders = {};
else
    runFolders = sort(unique(tmp));
end
end

function allDirs = collect_subdirs(rootDir)
rootDir = char(rootDir);
if ~isfolder(rootDir)
    allDirs = {};
    return;
end

% Robust recursive listing using MATLAB's path generation.
% This avoids cell/type edge cases in manual queue traversal.
allPath = genpath(rootDir);
allDirs = strsplit(allPath, pathsep);
allDirs = allDirs(~cellfun('isempty', allDirs));
allDirs = unique(allDirs, 'stable');
end

function overlayData = init_overlay_data(fields)
overlayData = struct();
for i = 1:numel(fields)
    overlayData.(fields{i}) = struct();
end
end

function overlayData = append_overlay_data(overlayData, fields, trLabel, thetaDeg, solid, shell)
groupTag = regexp(trLabel, 'tr\d+', 'match', 'once');
if isempty(groupTag)
    groupTag = 'misc';
end

for i = 1:numel(fields)
    fname = fields{i};
    if ~ismember(fname, solid.Properties.VariableNames) || ~ismember(fname, shell.Properties.VariableNames)
        continue;
    end

    if ~isfield(overlayData.(fname), groupTag)
        overlayData.(fname).(groupTag) = struct('thetaDeg', [], 'solid', {{}}, 'shell', {{}});
    end
    overlayData.(fname).(groupTag).thetaDeg(end+1) = thetaDeg;
    overlayData.(fname).(groupTag).solid{end+1} = double(solid.(fname));
    overlayData.(fname).(groupTag).shell{end+1} = double(shell.(fname));
end
end

function [trLabel, thetaDeg] = parse_run_tags_local(runName, runFolder)
trLabel = runName;
thetaDeg = NaN;

trToken = regexp(runName, 'tr\d+', 'match', 'once');
if ~isempty(trToken)
    trLabel = trToken;
end

angToken = regexp(runName, 'ang\d+', 'match', 'once');
if ~isempty(angToken)
    thetaDeg = str2double(angToken(4:end));
    return;
end

runNameLower = lower(runName);
if contains(runNameLower, 'dirx')
    thetaDeg = 0;
    return;
elseif contains(runNameLower, 'diry')
    thetaDeg = 90;
    return;
end

logPath = fullfile(runFolder, 'run_log.txt');
if ~isfile(logPath)
    return;
end

fid = fopen(logPath, 'r');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    lineLow = lower(strtrim(line));
    if startsWith(lineLow, 'line_direction:')
        if contains(lineLow, 'x')
            thetaDeg = 0;
        elseif contains(lineLow, 'y')
            thetaDeg = 90;
        end
        break;
    end
end
end

function [solid, shell] = load_run_pair_local(runFolder)
solid = load_spinodal_model(runFolder, 'solid');
shell = load_spinodal_model(runFolder, 'shell');

vars = solid.Properties.VariableNames;
for i = 1:numel(vars)
    v = solid.(vars{i});
    if ~isnumeric(v)
        solid.(vars{i}) = double(str2double(string(v)));
    else
        solid.(vars{i}) = double(v);
    end
end

vars = shell.Properties.VariableNames;
for i = 1:numel(vars)
    v = shell.(vars{i});
    if ~isnumeric(v)
        shell.(vars{i}) = double(str2double(string(v)));
    else
        shell.(vars{i}) = double(v);
    end
end

solid = sortrows(solid, 'Label');
shell = sortrows(shell, 'Label');

if ~isequal(solid.Label, shell.Label)
    error('run_all_analysis:LabelMismatch', ...
        'Label mismatch for %s (solid/shell nodes do not match).', runFolder);
end
end
