function nodalTable = read_abaqus_nodal_csv(csvPath)
%READ_ABAQUS_NODAL_CSV Read Abaqus nodal CSV and normalize Label column.
%
% nodalTable = READ_ABAQUS_NODAL_CSV(csvPath)
%
% Requirements:
% - CSV must contain a node label column (case-insensitive match to 'Label').
% - Label values are converted to integer doubles.
%
% Duplicate labels:
% - Kept deterministically using first occurrence.

    arguments
        csvPath (1, :) char
    end

    if ~isfile(csvPath)
        error('read_abaqus_nodal_csv:FileNotFound', 'File not found: %s', csvPath);
    end

    opts = detectImportOptions(csvPath);
    nodalTable = readtable(csvPath, opts);

    if isempty(nodalTable)
        error('read_abaqus_nodal_csv:EmptyCSV', 'CSV is empty: %s', csvPath);
    end

    varNames = nodalTable.Properties.VariableNames;
    labelIdx = find(strcmpi(varNames, 'Label'), 1);
    if isempty(labelIdx)
        error('read_abaqus_nodal_csv:MissingLabel', ...
            'CSV must contain a ''Label'' column (case-insensitive): %s', csvPath);
    end

    rawLabels = nodalTable.(varNames{labelIdx});
    labels = to_numeric_vector(rawLabels, varNames{labelIdx});
    if any(isnan(labels))
        error('read_abaqus_nodal_csv:BadLabelValues', ...
            'Label column contains non-numeric values: %s', csvPath);
    end
    labels = round(labels(:));

    % Place canonical Label column first.
    oldLabelName = varNames{labelIdx};
    nodalTable.(oldLabelName) = labels;
    if ~strcmp(oldLabelName, 'Label')
        nodalTable.Properties.VariableNames{labelIdx} = 'Label';
    end
    if labelIdx ~= 1
        nodalTable = movevars(nodalTable, 'Label', 'Before', 1);
    end

    % Remove duplicate labels deterministically (first row wins).
    [~, firstIdx] = unique(nodalTable.Label, 'stable');
    if numel(firstIdx) ~= height(nodalTable)
        nodalTable = nodalTable(firstIdx, :);
    end
end

function v = to_numeric_vector(x, nameForErrors)
    if isnumeric(x)
        v = double(x(:));
        return;
    end

    if iscellstr(x) || isstring(x) || iscell(x)
        v = str2double(string(x(:)));
        return;
    end

    error('read_abaqus_nodal_csv:UnsupportedType', ...
        'Unsupported type for column %s', nameForErrors);
end
