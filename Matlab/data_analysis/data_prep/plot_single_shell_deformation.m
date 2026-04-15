function [fig, profile] = plot_single_shell_deformation(csvPath)
%PLOT_SINGLE_SHELL_DEFORMATION Plot Z and U3 against Y for a shell CSV.
%   plot_single_shell_deformation() loads the default shell CSV and plots
%   Y on the x-axis with Z and U3 shown as separate curves.
%
%   [fig, profile] = plot_single_shell_deformation(csvPath) uses a custom
%   CSV path and returns the figure handle plus the processed profile data.

if nargin < 1 || isempty(csvPath)
    csvPath = default_csv_path();
end
csvPath = resolve_csv_path(char(csvPath));

if ~isfile(csvPath)
    error('plot_single_shell_deformation:MissingCsv', ...
        'CSV file not found: %s', csvPath);
end

tbl = readtable(csvPath);
requiredVars = {'Y', 'Z', 'U3'};
missingVars = requiredVars(~ismember(requiredVars, tbl.Properties.VariableNames));
if ~isempty(missingVars)
    error('plot_single_shell_deformation:MissingColumns', ...
        'CSV is missing required columns: %s', strjoin(missingVars, ', '));
end

Y = double(tbl.Y);
Z = double(tbl.Z);
U3 = double(tbl.U3);

valid = isfinite(Y) & isfinite(Z) & isfinite(U3);
if nnz(valid) < 2
    error('plot_single_shell_deformation:InsufficientData', ...
        'Need at least two valid rows with finite Y, Z, and U3 values.');
end

Y = Y(valid);
Z = Z(valid);
U3 = U3(valid);

[Y, order] = sort(Y);
Z = Z(order);
U3 = U3(order);

fig = figure('Color', 'w');
ax = axes(fig); %#ok<LAXES>
plot(ax, Y, Z, '-k', 'LineWidth', 1.5, 'DisplayName', 'Z');
hold(ax, 'on');
plot(ax, Y, U3, '-b', 'LineWidth', 1.5, 'DisplayName', 'U3');
grid(ax, 'on');
xlabel(ax, 'Y');
ylabel(ax, 'Z / U3');
xlim(ax, [0, 0.04]);
title(ax, 'Single shell Z and U3 profile', 'FontWeight', 'bold');
legend(ax, 'Location', 'best');

[parentDir, csvName, csvExt] = fileparts(csvPath);
[~, runName] = fileparts(parentDir);
subtitle(ax, sprintf('%s/%s%s', runName, csvName, csvExt), 'Interpreter', 'none');

profile = struct( ...
    'csvPath', csvPath, ...
    'Y', Y, ...
    'Z', Z, ...
    'U3', U3);

end

function csvPath = default_csv_path()
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');
csvPath = fullfile(repoRoot, ...
    'Matlab', 'GNN', 'data', 'dataset', 'samples', ...
    'sheetCone_tr50_ang030_lamellar_N128_1x1', ...
    'midpoint_results_shell.csv');
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
