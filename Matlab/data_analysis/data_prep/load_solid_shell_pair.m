function [solid, shell] = load_solid_shell_pair(trLabel, thetaDeg)
%LOAD_SOLID_SHELL_PAIR Load and align solid/shell data for a given ratio/angle.
%   [solid, shell] = load_solid_shell_pair(trLabel, thetaDeg) builds the
%   sheetCone folder path from the ratio (e.g. 'tr50') and lamellar angle
%   in degrees (e.g. 30), reads both CSV files, sorts them by Label, and
%   asserts that node identity and coordinates match within a tight tolerance.

arguments
    trLabel (1, :) char
    thetaDeg (1, 1) double
end

thetaTag = sprintf('ang%03d', round(thetaDeg));
modelFolder = sprintf('sheetCone_%s_%s_lamellar_N64_2x3', trLabel, thetaTag);

% Resolve paths relative to this function to avoid reliance on cwd.
thisDir = fileparts(mfilename('fullpath'));
% Raw data live under Matlab/results/..., which is two levels above this
% data_prep folder.
dataRoot = fullfile(thisDir, '..', '..', 'results', 'sheets', 'lamellar');
sheetFolder = fullfile(dataRoot, trLabel, modelFolder);
if ~isfolder(sheetFolder)
    error('load_solid_shell_pair:MissingFolder', ...
        'Could not find folder for %s at %s', modelFolder, sheetFolder);
end

solid = load_spinodal_model(sheetFolder, 'solid');
shell = load_spinodal_model(sheetFolder, 'shell');

% Force every column to numeric doubles to avoid downstream cell issues.
solid = force_all_numeric(solid);
shell = force_all_numeric(shell);

% Validate key columns
coordFields = {'X', 'Y', 'Z'};
allFields = [{'Label'}, coordFields];
for i = 1:numel(allFields)
    name = allFields{i};
    if ~isnumeric(solid.(name)) || ~isnumeric(shell.(name))
        error('load_solid_shell_pair:NonNumeric', ...
            'Column %s is not numeric after conversion (solid: %s, shell: %s)', ...
            name, class(solid.(name)), class(shell.(name)));
    end
end

solid = sortrows(solid, 'Label');
shell = sortrows(shell, 'Label');

if ~isequal(solid.Label, shell.Label)
    labelDiff = setxor(solid.Label, shell.Label);
    error('load_solid_shell_pair:LabelMismatch', ...
        'Label mismatch for %s: non-shared labels %s', ...
        modelFolder, mat2str(labelDiff(:).'));
end

coordFields = {'X', 'Y', 'Z'};
coordDelta = abs(table2array(solid(:, coordFields)) - table2array(shell(:, coordFields)));
maxCoordDelta = max(coordDelta, [], 'all');
tol = 1e-8;
if maxCoordDelta > tol
    error('load_solid_shell_pair:CoordinateMismatch', ...
        'Coordinate mismatch for %s: max |dXYZ| = %.3g exceeds tolerance %.1e', ...
        modelFolder, maxCoordDelta, tol);
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

end
