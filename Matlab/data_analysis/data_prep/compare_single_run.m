function [metrics, errorTable, outputDir] = compare_single_run(runFolder, opts)
%COMPARE_SINGLE_RUN Single-run comparison: overlays + curvature checks.
%   [metrics, errorTable, outputDir] = compare_single_run(runFolder)
%   loads FEA (solid) and FEA_shell (shell), generates overlay plots, and
%   runs all curvature comparisons.

arguments
    runFolder (1, :) char
    opts.OutputDir (1, :) char = ''
    opts.DoPlots (1, 1) logical = true
    opts.WriteFiles (1, 1) logical = true
end

runFolder = resolve_run_folder(char(runFolder));
if ~isfolder(runFolder)
    error('compare_single_run:MissingFolder', 'Run folder not found: %s', runFolder);
end

[~, runName] = fileparts(runFolder);

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');
resultsRoot = fullfile(scriptDir, '..', '..', 'results', 'analysis', 'single_runs');
if isempty(opts.OutputDir)
    outputDir = fullfile(resultsRoot, runName);
else
    outputDir = resolve_path(char(opts.OutputDir), repoRoot);
end

[solid, shell] = load_run_pair(runFolder);
[metrics, errorTable] = compute_basic_metrics(solid, shell);
[curvatureCheck, curvatureTable] = check_curvature_single_run(runFolder, solid, shell);

if opts.WriteFiles || opts.DoPlots
    ensure_dir(outputDir);
end

if opts.WriteFiles
    writetable(curvatureTable, fullfile(outputDir, 'curvature_check.csv'));
    save(fullfile(outputDir, 'curvature_check.mat'), 'curvatureCheck');
end

if opts.DoPlots
    [trLabel, thetaDeg] = parse_run_tags(runName, runFolder);
    plot_scatter_fields(solid, shell, trLabel, thetaDeg, outputDir);
    plot_curvature(solid, shell, trLabel, thetaDeg, outputDir, ...
        'CurvatureCheck', curvatureCheck);
    plot_curvature_offset_removed(solid, shell, trLabel, thetaDeg, outputDir, ...
        'CurvatureCheck', curvatureCheck);
end

if ~isnan(curvatureCheck.theory.kappa)
    fprintf(['[%s] Curvature check: theory kappa=%.6g 1/m | solid=%.6g 1/m ' ...
        '(rel=%.3g) | shell=%.6g 1/m (rel=%.3g)\n'], ...
        runName, curvatureCheck.theory.kappa, curvatureCheck.solid.kappaCircle, ...
        curvatureCheck.solid.RelErrTheory, curvatureCheck.shell.kappaCircle, ...
        curvatureCheck.shell.RelErrTheory);
else
    fprintf('[%s] Curvature check: run_log inputs not found (theory skipped).\n', runName);
end

fprintf('[%s] Comparison complete. Output: %s\n', runName, outputDir);

end

function runFolderOut = resolve_run_folder(runFolderIn)
runFolderOut = runFolderIn;
if isfolder(runFolderOut)
    return;
end

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');
alt = fullfile(repoRoot, runFolderIn);
if isfolder(alt)
    runFolderOut = alt;
end
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

function [trLabel, thetaDeg] = parse_run_tags(runName, runFolder)
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

% Laminate naming convention: dirx -> 0 deg, diry -> 90 deg.
runNameLower = lower(runName);
if contains(runNameLower, 'dirx')
    thetaDeg = 0;
    return;
elseif contains(runNameLower, 'diry')
    thetaDeg = 90;
    return;
end

% Fallback to run_log line_direction when available.
logPath = fullfile(runFolder, 'run_log.txt');
if isfile(logPath)
    fid = fopen(logPath, 'r');
    if fid >= 0
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
end
end

function [solid, shell] = load_run_pair(runFolder)
solid = load_spinodal_model(runFolder, 'solid');
shell = load_spinodal_model(runFolder, 'shell');

solid = force_all_numeric(solid);
shell = force_all_numeric(shell);

solid = sortrows(solid, 'Label');
shell = sortrows(shell, 'Label');

if ~isequal(solid.Label, shell.Label)
    labelDiff = setxor(solid.Label, shell.Label);
    error('compare_single_run:LabelMismatch', ...
        'Label mismatch for %s: non-shared labels %s', ...
        runFolder, mat2str(labelDiff(:).'));
end

coordFields = {'X', 'Y', 'Z'};
coordDelta = abs(table2array(solid(:, coordFields)) - table2array(shell(:, coordFields)));
maxCoordDelta = max(coordDelta, [], 'all');
tol = 1e-8;
if maxCoordDelta > tol
    error('compare_single_run:CoordinateMismatch', ...
        'Coordinate mismatch for %s: max |dXYZ| = %.3g exceeds tolerance %.1e', ...
        runFolder, maxCoordDelta, tol);
end
end

function tblOut = force_all_numeric(tblIn)
vars = tblIn.Properties.VariableNames;
tblOut = tblIn;
for i = 1:numel(vars)
    name = vars{i};
    v = tblOut.(name);
    if ~isnumeric(v)
        v = str2double(string(v));
    end
    tblOut.(name) = double(v);
end
end

function [metrics, errorTable] = compute_basic_metrics(solid, shell)
fields = {'U1','U2','U3','S_11','S_22','S_33','S_12','S_13','S_23','S_Mises'};

metrics = struct();
errorTable = table(double(solid.Label), double(solid.X), double(solid.Y), double(solid.Z), ...
    'VariableNames', {'Label', 'X', 'Y', 'Z'});

for i = 1:numel(fields)
    name = fields{i};
    if ~ismember(name, solid.Properties.VariableNames) || ~ismember(name, shell.Properties.VariableNames)
        continue;
    end
    s = double(solid.(name));
    sh = double(shell.(name));
    d = sh - s;
    mae = mean(abs(d), 'omitnan');
    rVal = NaN;
    if std(s) > 0 && std(sh) > 0
        cc = corrcoef(s, sh);
        if size(cc, 1) >= 2
            rVal = cc(1, 2);
        end
    end
    metrics.(name) = struct('MAE', mae, 'R', rVal);
    errorTable.(sprintf('d_%s', name)) = d;
end
end
