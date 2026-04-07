function tbl = load_spinodal_model(sheetFolder, modelType)
%LOAD_SPINODAL_MODEL Load a single spinodal model table from the FEA output.
%   tbl = load_spinodal_model(sheetFolder, modelType) reads the lone CSV file
%   from either the FEA (solid) or FEA_shell (shell) subfolder beneath
%   sheetFolder. The table is returned with a fixed variable order so that
%   downstream comparisons are deterministic.

arguments
    sheetFolder (1, :) char
    modelType (1, :) char
end

modelType = validatestring(modelType, {'solid', 'shell'});

subDir = 'FEA';
if strcmpi(modelType, 'shell')
    subDir = 'FEA_shell';
end

dataDir = fullfile(sheetFolder, subDir);
if ~isfolder(dataDir)
    error('load_spinodal_model:MissingFolder', ...
        'Folder not found for %s model: %s', modelType, dataDir);
end

csvFiles = dir(fullfile(dataDir, '*.csv'));
if isempty(csvFiles)
    error('load_spinodal_model:MissingCSV', ...
        'No CSV file found in %s', dataDir);
elseif numel(csvFiles) > 1
    error('load_spinodal_model:MultipleCSV', ...
        'Expected exactly one CSV in %s, found %d: %s', dataDir, numel(csvFiles), ...
        strjoin({csvFiles.name}, ', '));
end

filePath = fullfile(dataDir, csvFiles(1).name);

% Force numeric import for all expected columns to avoid cell arrays.
expectedVars = {'Label','X','Y','Z','U1','U2','U3','S_11','S_22','S_33','S_12','S_13','S_23','S_Mises'};

opts = detectImportOptions(filePath, 'Delimiter', ',');
tbl = readtable(filePath, opts);

% Normalize common stress column variants (e.g., S11 -> S_11).
renameMap = containers.Map( ...
    {'S11','S22','S33','S12','S13','S23','SMises'}, ...
    {'S_11','S_22','S_33','S_12','S_13','S_23','S_Mises'});
vars = tbl.Properties.VariableNames;
for iVar = 1:numel(vars)
    name = vars{iVar};
    if isKey(renameMap, name)
        vars{iVar} = renameMap(name);
    end
end
tbl.Properties.VariableNames = vars;

missingVars = setdiff(expectedVars, tbl.Properties.VariableNames);
if ~isempty(missingVars)
    error('load_spinodal_model:MissingColumns', ...
        'File %s is missing required columns: %s', filePath, strjoin(missingVars, ', '));
end

tbl = tbl(:, expectedVars);
for iVar = 1:numel(expectedVars)
    name = expectedVars{iVar};
    v = tbl.(name);
    if iscell(v)
        v = str2double(string(v));
    elseif isstring(v)
        v = str2double(v);
    end
    tbl.(name) = double(v);
end

end
