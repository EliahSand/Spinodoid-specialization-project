function boundary = detect_boundary_nodes(elementLocalNodes, numNodes)
%DETECT_BOUNDARY_NODES Find boundary nodes from shell element connectivity.
%
% boundary = DETECT_BOUNDARY_NODES(elementLocalNodes, numNodes)
%
% Boundary criterion:
% - An edge is boundary if it is used by exactly one element.
% - Boundary nodes are nodes belonging to at least one boundary edge.
%
% Output fields:
%   all_edges        - all extracted undirected edges (with duplicates)
%   unique_edges     - unique undirected edges
%   edge_use_count   - multiplicity for each unique edge
%   boundary_edges   - unique edges with use-count == 1
%   boundary_mask    - numNodes x 1 logical

    arguments
        elementLocalNodes
        numNodes (1, 1) double {mustBeNonnegative, mustBeInteger}
    end

    [allEdges, ~] = extract_shell_edges(elementLocalNodes);
    if isempty(allEdges)
        boundary = struct();
        boundary.all_edges = zeros(0, 2);
        boundary.unique_edges = zeros(0, 2);
        boundary.edge_use_count = zeros(0, 1);
        boundary.boundary_edges = zeros(0, 2);
        boundary.boundary_mask = false(numNodes, 1);
        return;
    end

    [uniqueEdges, ~, edgeClass] = unique(allEdges, 'rows', 'stable');
    edgeUseCount = accumarray(edgeClass, 1);
    boundaryEdges = uniqueEdges(edgeUseCount == 1, :);

    boundaryMask = false(numNodes, 1);
    if ~isempty(boundaryEdges)
        boundaryMask(unique(boundaryEdges(:))) = true;
    end

    boundary = struct();
    boundary.all_edges = allEdges;
    boundary.unique_edges = uniqueEdges;
    boundary.edge_use_count = edgeUseCount;
    boundary.boundary_edges = boundaryEdges;
    boundary.boundary_mask = boundaryMask;
end

