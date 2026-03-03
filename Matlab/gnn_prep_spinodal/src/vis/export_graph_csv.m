function paths = export_graph_csv(graphData, outDir, prefix)
%EXPORT_GRAPH_CSV Export graph nodes/edges to CSV (+ MAT copy).
%
% paths = EXPORT_GRAPH_CSV(graphData, outDir, prefix)
%
% Outputs:
%   <prefix>_nodes.csv
%   <prefix>_edges.csv
%   <prefix>.mat

    arguments
        graphData struct
        outDir (1, :) char
        prefix (1, :) char = 'graph'
    end

    if ~isfolder(outDir)
        mkdir(outDir);
    end

    n = numel(graphData.node_labels);
    coords = graphData.node_coords;
    boundaryMask = false(n, 1);
    if isfield(graphData, 'boundary_mask')
        b = logical(graphData.boundary_mask(:));
        if numel(b) == n
            boundaryMask = b;
        end
    end

    nodeTable = table( ...
        (1:n).', ...
        graphData.node_labels(:), ...
        coords(:, 1), coords(:, 2), coords(:, 3), ...
        boundaryMask, ...
        'VariableNames', {'NodeIndex', 'Label', 'X', 'Y', 'Z', 'IsBoundary'});

    if isfield(graphData, 'node_data_table') && istable(graphData.node_data_table) ...
            && height(graphData.node_data_table) == n
        attrTable = graphData.node_data_table;
        if any(strcmp(attrTable.Properties.VariableNames, 'Label'))
            attrTable(:, 'Label') = [];
        end

        % Avoid variable-name collisions (e.g., CSV has X/Y/Z while nodeTable already has X/Y/Z).
        usedNames = nodeTable.Properties.VariableNames;
        attrNames = attrTable.Properties.VariableNames;
        newAttrNames = attrNames;
        for i = 1:numel(attrNames)
            base = attrNames{i};
            candidate = base;
            if any(strcmp(candidate, usedNames)) || any(strcmp(candidate, newAttrNames(1:i-1)))
                candidate = ['Attr_' base];
                k = 1;
                while any(strcmp(candidate, usedNames)) || any(strcmp(candidate, newAttrNames(1:i-1)))
                    candidate = sprintf('Attr_%s_%d', base, k);
                    k = k + 1;
                end
            end
            newAttrNames{i} = candidate;
            usedNames{end + 1} = candidate; %#ok<AGROW>
        end
        attrTable.Properties.VariableNames = newAttrNames;

        nodeTable = [nodeTable, attrTable];
    end

    edges = graphData.edges_local;
    if isempty(edges)
        edgeTable = table('Size', [0, 4], ...
            'VariableTypes', {'double', 'double', 'double', 'double'}, ...
            'VariableNames', {'SourceIndex', 'TargetIndex', 'SourceLabel', 'TargetLabel'});
    else
        src = edges(:, 1);
        dst = edges(:, 2);
        edgeTable = table( ...
            src, dst, ...
            graphData.node_labels(src), ...
            graphData.node_labels(dst), ...
            'VariableNames', {'SourceIndex', 'TargetIndex', 'SourceLabel', 'TargetLabel'});
    end

    nodesPath = fullfile(outDir, sprintf('%s_nodes.csv', prefix));
    edgesPath = fullfile(outDir, sprintf('%s_edges.csv', prefix));
    matPath = fullfile(outDir, sprintf('%s.mat', prefix));

    writetable(nodeTable, nodesPath);
    writetable(edgeTable, edgesPath);
    save(matPath, 'graphData', '-v7.3');

    paths = struct();
    paths.nodes_csv_path = nodesPath;
    paths.edges_csv_path = edgesPath;
    paths.mat_path = matPath;
end
