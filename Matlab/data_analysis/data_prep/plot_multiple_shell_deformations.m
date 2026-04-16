function [fig, profiles] = plot_multiple_shell_deformations(csvPaths)
%PLOT_MULTIPLE_SHELL_DEFORMATIONS Overlay multiple shell midpoint CSV runs.
%   plot_multiple_shell_deformations() loads a hardcoded list of shell CSV
%   files and overlays their U3(Y) curves in a single figure.
%
%   [fig, profiles] = plot_multiple_shell_deformations(csvPaths) uses a
%   custom cell array / string array of CSV paths and returns the figure
%   handle plus the processed profile data.

if nargin < 1 || isempty(csvPaths)
    csvPaths = default_csv_paths();
end

if ischar(csvPaths) || isstring(csvPaths)
    csvPaths = cellstr(csvPaths);
end

if ~iscell(csvPaths) || isempty(csvPaths)
    error('plot_multiple_shell_deformations:BadInput', ...
        'csvPaths must be a non-empty cell array, string array, or char path.');
end

nRuns = numel(csvPaths);
profiles = repmat(struct( ...
    'csvPath', '', ...
    'label', '', ...
    'Y', [], ...
    'U3', []), 1, nRuns);

fig = figure('Color', 'w');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
grid(ax, 'on');
xlabel(ax, 'Y');
ylabel(ax, 'U3');
xlim(ax, [0, 0.04]);
title(ax, sprintf('Multiple shell U3 overlays (%d runs)', nRuns), ...
    'FontWeight', 'bold');
yline(ax, 0, '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');

colors = lines(nRuns);
for i = 1:nRuns
    profile = load_shell_profile(csvPaths{i});
    profiles(i) = profile;

    plot(ax, profile.Y, profile.U3, 'LineWidth', 1.5, ...
        'Color', colors(i, :), 'DisplayName', profile.label);
end

legend(ax, 'Location', 'best', 'Interpreter', 'none');

end

function csvPaths = default_csv_paths()
csvPaths = { ...
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang000_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang010_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang020_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang030_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang040_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang050_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang060_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang070_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang080_lamellar_N128_1x1/midpoint_results_shell.csv'
    'Matlab/GNN/data/dataset/samples/sheetCone_tr50_ang090_lamellar_N128_1x1/midpoint_results_shell.csv'
    };
end

function profile = load_shell_profile(csvPath)
csvPath = resolve_csv_path(char(csvPath));
if ~isfile(csvPath)
    error('plot_multiple_shell_deformations:MissingCsv', ...
        'CSV file not found: %s', csvPath);
end

tbl = readtable(csvPath);
requiredVars = {'Y', 'U3'};
missingVars = requiredVars(~ismember(requiredVars, tbl.Properties.VariableNames));
if ~isempty(missingVars)
    error('plot_multiple_shell_deformations:MissingColumns', ...
        'CSV %s is missing required columns: %s', ...
        csvPath, strjoin(missingVars, ', '));
end

Y = double(tbl.Y);
U3 = double(tbl.U3);

valid = isfinite(Y) & isfinite(U3);
if nnz(valid) < 2
    error('plot_multiple_shell_deformations:InsufficientData', ...
        'CSV %s needs at least two valid rows with finite Y and U3 values.', ...
        csvPath);
end

Y = Y(valid);
U3 = U3(valid);

[Y, order] = sort(Y);
U3 = U3(order);

[parentDir, ~, ~] = fileparts(csvPath);
[~, runName] = fileparts(parentDir);

profile = struct( ...
    'csvPath', csvPath, ...
    'label', runName, ...
    'Y', Y, ...
    'U3', U3);
end

function pOut = resolve_csv_path(pIn)
pOut = pIn;
if isempty(pIn) || is_absolute_path(pIn)
    return;
end

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');
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
