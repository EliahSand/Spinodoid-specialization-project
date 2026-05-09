function [matGraph, debugData] = extract_mat_spline_graph_gnn(fullGraph, varargin)
%EXTRACT_MAT_SPLINE_GRAPH_GNN Build MAT/spline control graph from full mesh graph.
%
% The dense shell mesh nodes are rasterized back to their native XY grid.
% The MAT/spline graph is then built from the occupied mask and returned as
% a compact control-point graph for GNN schema v4.

    required = {'node_coords', 'node_labels', 'edges_local', 'boundary_mask'};
    for i = 1:numel(required)
        if ~isfield(fullGraph, required{i})
            error('extract_mat_spline_graph_gnn:MissingField', ...
                'fullGraph is missing required field "%s".', required{i});
        end
    end

    coords = fullGraph.node_coords;
    if size(coords, 2) < 2
        error('extract_mat_spline_graph_gnn:BadCoords', ...
            'fullGraph.node_coords must have at least two columns [X,Y].');
    end

    xy = double(coords(:, 1:2));
    [xVals, xBin] = cluster_axis_values(xy(:, 1));
    [yVals, yBin] = cluster_axis_values(xy(:, 2));

    nodeAt = zeros(numel(yVals), numel(xVals));
    for i = 1:size(xy, 1)
        r = yBin(i);
        c = xBin(i);
        if nodeAt(r, c) == 0
            nodeAt(r, c) = i;
        end
    end
    occMask = nodeAt > 0;

    [matGraph, debugData] = mat_spline_graph_from_mask(occMask, xVals, yVals, varargin{:});
    debugData.node_at_grid = nodeAt;
    debugData.x_grid = xVals;
    debugData.y_grid = yVals;

    matGraph.parent_kind = 'full_reference_graph';
    matGraph.source_num_nodes = fullGraph.num_nodes;
    matGraph.source_num_edges = size(fullGraph.edges_local, 1);
    if isfield(fullGraph, 'source_inp_path')
        matGraph.source_inp_path = fullGraph.source_inp_path;
    end
end

function [valsOut, binIdx] = cluster_axis_values(vals)
    vals = double(vals(:));
    n = numel(vals);
    if n == 0
        valsOut = zeros(0, 1);
        binIdx = zeros(0, 1);
        return;
    end

    [sortedVals, order] = sort(vals);
    diffs = diff(sortedVals);
    posDiffs = diffs(diffs > 0);
    if isempty(posDiffs)
        tol = max(1e-12 * max(abs(sortedVals)), 1e-15);
    else
        tol = max(min(posDiffs) * 0.25, 1e-12 * max([range(sortedVals), 1]));
    end

    groupSorted = ones(n, 1);
    g = 1;
    for i = 2:n
        if abs(sortedVals(i) - sortedVals(i - 1)) > tol
            g = g + 1;
        end
        groupSorted(i) = g;
    end

    valsOut = accumarray(groupSorted, sortedVals, [], @mean);
    binIdx = zeros(n, 1);
    binIdx(order) = groupSorted;
end
