function graphOut = attach_nodal_data(graphIn, nodalTable, varargin)
%ATTACH_NODAL_DATA Map nodal CSV fields to graph nodes by Abaqus label.
%
% graphOut = ATTACH_NODAL_DATA(graphIn, nodalTable, ...)
%
% Name-Value options:
%   'Strict'    : if true, error when a graph node has no CSV row (default true)
%   'FillValue' : value used for missing numeric fields when Strict=false (default NaN)
%
% Added fields in graphOut:
%   node_data_table      - table aligned with graph node ordering
%   node_attr_names      - attribute names (excluding Label)
%   csv_mapping          - mapping metadata

    p = inputParser;
    p.addParameter('Strict', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('FillValue', NaN, @isnumeric);
    p.parse(varargin{:});
    opts = p.Results;

    if ~istable(nodalTable)
        error('attach_nodal_data:BadCSVTable', 'nodalTable must be a MATLAB table.');
    end

    varNames = nodalTable.Properties.VariableNames;
    labelIdx = find(strcmpi(varNames, 'Label'), 1);
    if isempty(labelIdx)
        error('attach_nodal_data:MissingLabel', 'nodalTable must contain a Label column.');
    end

    csvLabels = round(double(nodalTable{:, labelIdx}));
    if any(isnan(csvLabels))
        error('attach_nodal_data:BadLabels', 'Label column contains non-numeric values.');
    end

    % Remove duplicate labels deterministically (first row wins).
    [~, firstRows] = unique(csvLabels, 'stable');
    if numel(firstRows) ~= numel(csvLabels)
        nodalTable = nodalTable(firstRows, :);
        csvLabels = csvLabels(firstRows);
    end

    graphLabels = graphIn.node_labels(:);
    [isMapped, csvRowForNode] = ismember(graphLabels, csvLabels);
    if logical(opts.Strict) && any(~isMapped)
        missingCount = nnz(~isMapped);
        error('attach_nodal_data:MissingNodeRows', ...
            '%d graph nodes are missing from nodal CSV labels.', missingCount);
    end

    attrNames = nodalTable.Properties.VariableNames;
    attrNames(labelIdx) = [];

    nNodes = numel(graphLabels);
    aligned = table(graphLabels, 'VariableNames', {'Label'});
    for i = 1:numel(attrNames)
        name = attrNames{i};
        src = nodalTable.(name);
        aligned.(name) = align_column(src, csvRowForNode, isMapped, nNodes, opts.FillValue);
    end

    graphOut = graphIn;
    graphOut.node_data_table = aligned;
    graphOut.node_attr_names = attrNames;
    graphOut.csv_mapping = struct( ...
        'csv_row_for_node_index', csvRowForNode, ...
        'is_mapped', isMapped, ...
        'num_mapped', nnz(isMapped), ...
        'num_unmapped', nnz(~isMapped), ...
        'source_csv_labels', csvLabels);
end

function out = align_column(src, csvRowForNode, isMapped, nNodes, fillValue)
    if isnumeric(src) || islogical(src)
        out = repmat(fillValue, nNodes, 1);
        out(isMapped) = double(src(csvRowForNode(isMapped)));
        return;
    end

    if isstring(src)
        out = strings(nNodes, 1);
        out(isMapped) = src(csvRowForNode(isMapped));
        return;
    end

    if iscell(src)
        out = repmat({[]}, nNodes, 1);
        out(isMapped) = src(csvRowForNode(isMapped));
        return;
    end

    error('attach_nodal_data:UnsupportedColumnType', ...
        'Unsupported nodal attribute type: %s', class(src));
end

