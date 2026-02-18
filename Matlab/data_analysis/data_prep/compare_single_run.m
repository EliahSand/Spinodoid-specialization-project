function [metrics, errorTable, outputDir] = compare_single_run(runFolder, opts)
%COMPARE_SINGLE_RUN Compare SOLID vs FEA_shell for a single run folder.
%   [metrics, errorTable, outputDir] = compare_single_run(runFolder) loads
%   FEA (solid) and FEA_shell (shell) CSVs from runFolder, computes node-wise
%   errors, and writes analysis outputs under results/analysis/single_runs.
%
%   Optional name-value pairs:
%     OutputDir  - override output directory (default: results/analysis/single_runs/<runName>)
%     DoPlots    - enable diagnostic plots (default: true)
%     DoPCA      - enable PCA error analysis (default: false)
%     WriteFiles - write CSV/MAT outputs (default: true)

arguments
    runFolder (1, :) char
    opts.OutputDir (1, :) char = ''
    opts.DoPlots (1, 1) logical = true
    opts.DoPCA (1, 1) logical = false
    opts.WriteFiles (1, 1) logical = true
end

if ~isfolder(runFolder)
    error('compare_single_run:MissingFolder', 'Run folder not found: %s', runFolder);
end

runFolder = char(runFolder);
[~, runName] = fileparts(runFolder);

scriptDir = fileparts(mfilename('fullpath'));
resultsRoot = fullfile(scriptDir, '..', 'results', 'analysis', 'single_runs');
if isempty(opts.OutputDir)
    outputDir = fullfile(resultsRoot, runName);
else
    outputDir = opts.OutputDir;
end

[solid, shell] = load_run_pair(runFolder);
[metrics, errorTable] = compare_solid_shell(solid, shell);
[curvatureCheck, curvatureTable] = check_curvature_single_run(runFolder, solid, shell);

if opts.WriteFiles || opts.DoPlots
    ensure_dir(outputDir);
end

if opts.WriteFiles
    metricsTbl = metrics_to_table(metrics);
    writetable(metricsTbl, fullfile(outputDir, 'metrics_summary.csv'));
    writetable(errorTable, fullfile(outputDir, 'node_errors.csv'));
    writetable(curvatureTable, fullfile(outputDir, 'curvature_check.csv'));
    save(fullfile(outputDir, 'metrics.mat'), 'metrics');
    save(fullfile(outputDir, 'error_table.mat'), 'errorTable');
    save(fullfile(outputDir, 'curvature_check.mat'), 'curvatureCheck');
end

if opts.DoPlots
    [trLabel, thetaDeg] = parse_run_tags(runName, runFolder);
    plot_scatter_fields(solid, shell, metrics, trLabel, thetaDeg, outputDir);
    plot_curvature(solid, shell, trLabel, thetaDeg, outputDir, ...
        'CurvatureCheck', curvatureCheck);
    plot_curvature_offset_removed(solid, shell, trLabel, thetaDeg, outputDir, ...
        'CurvatureCheck', curvatureCheck);
    %plot_error_histograms(errorTable, metrics, trLabel, thetaDeg, outputDir);
    %plot_curvature_profiles(solid, shell, metrics, trLabel, thetaDeg, outputDir);
    if opts.DoPCA
        if license('test', 'statistics_toolbox') && exist('pca', 'file') == 2
            analyze_error_patterns(solid, shell, trLabel, thetaDeg, outputDir);
        else
            warning('compare_single_run:PCAUnavailable', ...
                'Skipping PCA: Statistics Toolbox (pca) not available.');
        end
    end
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
