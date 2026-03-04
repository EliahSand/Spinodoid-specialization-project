function [structGraph, skeletonGraph, debugData, annotatedFullGraph] = extract_structural_graph(fullGraph, varargin)
%EXTRACT_STRUCTURAL_GRAPH Build a simplified skeleton graph from a dense graph.
%
% [structGraph, skeletonGraph, debugData, annotatedFullGraph] = EXTRACT_STRUCTURAL_GRAPH(fullGraph, ...)
%
% Name-Value options:
%   'DetailLevel'       : [0,1], 0 = coarsest, 1 = keep full skeleton chains
%   'MinIslandNodes'    : minimum skeleton-node count kept for an isolated component
%
% The dense graph is first rasterized onto its native XY node grid, thinned to a
% one-node-thick skeleton with a standard Zhang-Suen topology-preserving thinning
% procedure, and then compressed by collapsing degree-2 chains. Reduced edges keep
% the full skeleton polyline and the originating dense-graph node ids.

    p = inputParser;
    p.addParameter('DetailLevel', 0.20, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('MinIslandNodes', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.parse(varargin{:});
    opts = p.Results;

    required = {'node_coords', 'node_labels', 'edges_local', 'boundary_mask'};
    for i = 1:numel(required)
        if ~isfield(fullGraph, required{i})
            error('extract_structural_graph:MissingField', ...
                'fullGraph is missing required field "%s".', required{i});
        end
    end

    coords = fullGraph.node_coords;
    if size(coords, 2) < 2
        error('extract_structural_graph:BadCoords', ...
            'fullGraph.node_coords must have at least two columns [X,Y].');
    end
    xy = coords(:, 1:2);
    nFull = size(xy, 1);

    [xVals, xBin] = cluster_axis_values(xy(:, 1));
    [yVals, yBin] = cluster_axis_values(xy(:, 2));
    nCols = numel(xVals);
    nRows = numel(yVals);

    nodeAt = zeros(nRows, nCols);
    for i = 1:nFull
        r = yBin(i);
        c = xBin(i);
        if nodeAt(r, c) == 0
            nodeAt(r, c) = i;
        end
    end
    occMask = nodeAt > 0;
    componentAt = zeros(nRows, nCols);
    if isfield(fullGraph, 'graph')
        fullCompId = conncomp(fullGraph.graph);
    else
        fullCompId = transpose(1:nFull);
    end
    for i = 1:nFull
        componentAt(yBin(i), xBin(i)) = fullCompId(i);
    end

    boundaryMask = false(nRows, nCols);
    fullBoundary = logical(fullGraph.boundary_mask(:));
    if numel(fullBoundary) ~= nFull
        error('extract_structural_graph:BadBoundaryMask', ...
            'boundary_mask length does not match node count.');
    end
    for i = find(fullBoundary(:)).'
        boundaryMask(yBin(i), xBin(i)) = true;
    end
    if ~any(boundaryMask(:))
        boundaryMask = find_perimeter_pixels(occMask);
    end

    boundaryStepDistance = multi_source_grid_distance(occMask, boundaryMask);
    [radiusMap, thicknessMap, spacingInfo] = compute_physical_radius_map(occMask, boundaryMask, xVals, yVals);
    annotatedFullGraph = attach_full_graph_thickness(fullGraph, xBin, yBin, radiusMap, thicknessMap);

    skelMask = false(size(occMask));
    for compId = 1:max(componentAt(:))
        compMask = componentAt == compId;
        if ~any(compMask(:))
            continue;
        end
        compSkel = zhang_suen_thinning(compMask);
        compSkel = enforce_component_representatives(compSkel, compMask, boundaryStepDistance, opts.MinIslandNodes);
        skelMask = skelMask | compSkel;
    end

    parentKind = get_default(annotatedFullGraph, 'kind', 'full_reference_graph');
    fullNodeData = table();
    if isfield(annotatedFullGraph, 'node_data_table') && istable(annotatedFullGraph.node_data_table) ...
            && height(annotatedFullGraph.node_data_table) == nFull
        fullNodeData = annotatedFullGraph.node_data_table;
    end

    [skeletonGraph, skelDebug] = build_skeleton_graph_from_mask( ...
        skelMask, componentAt, nodeAt, xy, annotatedFullGraph.node_labels(:), fullNodeData, boundaryStepDistance, parentKind, ...
        annotatedFullGraph.node_radius(:), annotatedFullGraph.node_thickness(:));

    structGraph = compress_skeleton_graph(skeletonGraph, opts.DetailLevel);
    structGraphBeforeCleanup = structGraph;
    [structGraph, cleanupDebug] = cleanup_structural_graph(structGraph, skeletonGraph);
    structGraph = attach_structural_edge_thickness(structGraph, annotatedFullGraph);

    debugData = struct();
    debugData.x_grid = xVals;
    debugData.y_grid = yVals;
    debugData.node_at_grid = nodeAt;
    debugData.occupancy_mask = occMask;
    debugData.boundary_mask = boundaryMask;
    debugData.boundary_step_distance = boundaryStepDistance;
    debugData.component_at_grid = componentAt;
    debugData.radius_map = radiusMap;
    debugData.thickness_map = thicknessMap;
    debugData.spacing_info = spacingInfo;
    debugData.skeleton_mask = skelMask;
    debugData.skeleton_local_at_grid = skelDebug.skelLocalAtGrid;
    debugData.structural_before_cleanup = structGraphBeforeCleanup;
    debugData.structural_cleanup = cleanupDebug;
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

function perim = find_perimeter_pixels(mask)
    [nRows, nCols] = size(mask);
    perim = false(nRows, nCols);
    [rows, cols] = find(mask);
    for k = 1:numel(rows)
        r = rows(k);
        c = cols(k);
        neigh = [r-1 c; r+1 c; r c-1; r c+1];
        for j = 1:size(neigh, 1)
            rr = neigh(j, 1);
            cc = neigh(j, 2);
            if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc)
                perim(r, c) = true;
                break;
            end
        end
    end
end

function distMap = multi_source_grid_distance(mask, seedMask)
    [nRows, nCols] = size(mask);
    distMap = inf(nRows, nCols);
    seedMask = logical(seedMask) & logical(mask);
    queue = zeros(nnz(mask), 2);
    qHead = 1;
    qTail = 0;

    [seedRows, seedCols] = find(seedMask);
    for k = 1:numel(seedRows)
        r = seedRows(k);
        c = seedCols(k);
        distMap(r, c) = 0;
        qTail = qTail + 1;
        queue(qTail, :) = [r, c];
    end

    if qTail == 0
        return;
    end

    neigh = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];
    while qHead <= qTail
        r = queue(qHead, 1);
        c = queue(qHead, 2);
        qHead = qHead + 1;
        base = distMap(r, c);

        for j = 1:size(neigh, 1)
            rr = r + neigh(j, 1);
            cc = c + neigh(j, 2);
            if rr < 1 || rr > nRows || cc < 1 || cc > nCols || ~mask(rr, cc)
                continue;
            end
            cand = base + 1;
            if cand < distMap(rr, cc)
                distMap(rr, cc) = cand;
                qTail = qTail + 1;
                queue(qTail, :) = [rr, cc];
            end
        end
    end
end

function [radiusMap, thicknessMap, spacingInfo] = compute_physical_radius_map(mask, boundaryMask, xVals, yVals)
    radiusMap = inf(size(mask));
    thicknessMap = inf(size(mask));

    dx = diff(sort(xVals(:)));
    dy = diff(sort(yVals(:)));
    dx = dx(dx > 0);
    dy = dy(dy > 0);
    spacingInfo = struct();
    spacingInfo.dx = median_or_nan(dx);
    spacingInfo.dy = median_or_nan(dy);
    spacingInfo.h = median_or_nan([dx(:); dy(:)]);

    seeds = find(mask & boundaryMask);
    if isempty(seeds)
        return;
    end

    occIdx = find(mask);
    occMask = false(size(mask));
    occMask(occIdx) = true;

    visited = false(size(mask));
    radiusMap(seeds) = 0;

    while true
        candidates = radiusMap;
        candidates(~occMask | visited) = inf;
        [bestDist, linearIdx] = min(candidates(:));
        if ~isfinite(bestDist)
            break;
        end

        visited(linearIdx) = true;
        [r, c] = ind2sub(size(mask), linearIdx);
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue;
                end
                rr = r + dr;
                cc = c + dc;
                if rr < 1 || rr > size(mask, 1) || cc < 1 || cc > size(mask, 2) || ~occMask(rr, cc)
                    continue;
                end
                step = hypot(xVals(cc) - xVals(c), yVals(rr) - yVals(r));
                cand = bestDist + step;
                if cand < radiusMap(rr, cc)
                    radiusMap(rr, cc) = cand;
                end
            end
        end
    end

    thicknessMap(mask) = 2 * radiusMap(mask);
end

function out = median_or_nan(vals)
    vals = vals(isfinite(vals));
    if isempty(vals)
        out = NaN;
    else
        out = median(vals);
    end
end

function graphOut = attach_full_graph_thickness(graphIn, xBin, yBin, radiusMap, thicknessMap)
    graphOut = graphIn;
    n = numel(graphIn.node_labels);
    radius = zeros(n, 1);
    thickness = zeros(n, 1);
    for i = 1:n
        radius(i) = radiusMap(yBin(i), xBin(i));
        thickness(i) = thicknessMap(yBin(i), xBin(i));
    end

    graphOut.node_radius = radius(:);
    graphOut.node_thickness = thickness(:);

    if ~isfield(graphOut, 'node_data_table') || ~istable(graphOut.node_data_table) || height(graphOut.node_data_table) ~= n
        graphOut.node_data_table = table(graphOut.node_labels(:), 'VariableNames', {'Label'});
    end
    graphOut.node_data_table = upsert_table_column(graphOut.node_data_table, 'Radius', radius(:));
    graphOut.node_data_table = upsert_table_column(graphOut.node_data_table, 'Thickness', thickness(:));
end

function skel = zhang_suen_thinning(mask)
    work = false(size(mask) + 2);
    work(2:(end - 1), 2:(end - 1)) = logical(mask);
    changed = true;
    while changed
        changed = false;

        toDelete = thinning_iteration(work, 1);
        if any(toDelete(:))
            work(toDelete) = false;
            changed = true;
        end

        toDelete = thinning_iteration(work, 2);
        if any(toDelete(:))
            work(toDelete) = false;
            changed = true;
        end
    end
    skel = work(2:(end - 1), 2:(end - 1));
end

function toDelete = thinning_iteration(skel, iterId)
    [nRows, nCols] = size(skel);
    toDelete = false(nRows, nCols);

    for r = 2:(nRows - 1)
        for c = 2:(nCols - 1)
            if ~skel(r, c)
                continue;
            end

            p2 = skel(r - 1, c);
            p3 = skel(r - 1, c + 1);
            p4 = skel(r, c + 1);
            p5 = skel(r + 1, c + 1);
            p6 = skel(r + 1, c);
            p7 = skel(r + 1, c - 1);
            p8 = skel(r, c - 1);
            p9 = skel(r - 1, c - 1);
            neigh = [p2 p3 p4 p5 p6 p7 p8 p9];

            n = sum(neigh);
            if n < 2 || n > 6
                continue;
            end

            transitions = sum(~neigh & [neigh(2:end), neigh(1)]);
            if transitions ~= 1
                continue;
            end

            if iterId == 1
                condA = ~(p2 && p4 && p6);
                condB = ~(p4 && p6 && p8);
            else
                condA = ~(p2 && p4 && p8);
                condB = ~(p2 && p6 && p8);
            end

            if condA && condB
                toDelete(r, c) = true;
            end
        end
    end
end

function skelOut = enforce_component_representatives(skelMask, occMask, distMap, minIslandNodes)
    skelOut = skelMask;

    occComponents = connected_components_binary(occMask, 8);
    for compId = 1:max(occComponents(:))
        compMask = occComponents == compId;
        if ~any(compMask(:))
            continue;
        end

        skelCount = nnz(skelOut & compMask);
        if skelCount >= minIslandNodes
            continue;
        end

        compIdx = find(compMask);
        compDist = distMap(compIdx);
        if all(isinf(compDist))
            [~, pickRel] = min(compIdx);
        else
            [~, pickRel] = max(compDist);
        end
        skelOut(compIdx(pickRel)) = true;
    end
end

function labels = connected_components_binary(mask, conn)
    if nargin < 2
        conn = 8;
    end
    [nRows, nCols] = size(mask);
    labels = zeros(nRows, nCols);
    current = 0;
    queue = zeros(nnz(mask), 2);
    if conn == 4
        neigh = [-1 0; 1 0; 0 -1; 0 1];
    else
        neigh = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];
    end

    [rows, cols] = find(mask);
    for seedIdx = 1:numel(rows)
        r0 = rows(seedIdx);
        c0 = cols(seedIdx);
        if labels(r0, c0) ~= 0
            continue;
        end

        current = current + 1;
        qHead = 1;
        qTail = 1;
        queue(1, :) = [r0, c0];
        labels(r0, c0) = current;

        while qHead <= qTail
            r = queue(qHead, 1);
            c = queue(qHead, 2);
            qHead = qHead + 1;

            for j = 1:size(neigh, 1)
                rr = r + neigh(j, 1);
                cc = c + neigh(j, 2);
                if rr < 1 || rr > nRows || cc < 1 || cc > nCols
                    continue;
                end
                if ~mask(rr, cc) || labels(rr, cc) ~= 0
                    continue;
                end
                qTail = qTail + 1;
                queue(qTail, :) = [rr, cc];
                labels(rr, cc) = current;
            end
        end
    end
end

function [skelGraph, debugData] = build_skeleton_graph_from_mask(skelMask, componentAt, nodeAt, xyFull, fullLabels, fullNodeData, distMap, parentKind, fullRadius, fullThickness)
    skelLinear = find(skelMask);
    nSkel = numel(skelLinear);
    if nSkel == 0
        error('extract_structural_graph:EmptySkeleton', ...
            'Skeleton extraction produced no nodes.');
    end

    fullIdx = nodeAt(skelLinear);
    coords = xyFull(fullIdx, :);

    skelLocalAt = zeros(size(skelMask));
    skelLocalAt(skelLinear) = 1:nSkel;

    dirs = [0 1; 1 -1; 1 0; 1 1];
    edges = zeros(0, 2);
    lengths = zeros(0, 1);
    for i = 1:nSkel
        [r, c] = ind2sub(size(skelMask), skelLinear(i));
        for d = 1:size(dirs, 1)
            rr = r + dirs(d, 1);
            cc = c + dirs(d, 2);
            if rr < 1 || rr > size(skelMask, 1) || cc < 1 || cc > size(skelMask, 2)
                continue;
            end
            thisComp = componentAt(r, c);
            if thisComp == 0 || componentAt(rr, cc) ~= thisComp
                continue;
            end
            if abs(dirs(d, 1)) == 1 && abs(dirs(d, 2)) == 1
                % Prevent corner-touch shortcuts across distinct spinodal regions.
                if ~supports_diagonal_link(componentAt, r, c, rr, cc, thisComp)
                    continue;
                end
            end
            j = skelLocalAt(rr, cc);
            if j == 0
                continue;
            end
            edges(end + 1, :) = [i, j]; %#ok<AGROW>
            lengths(end + 1, 1) = norm(coords(i, :) - coords(j, :)); %#ok<AGROW>
        end
    end

    if isempty(edges)
        G = graph([], [], [], nSkel);
    else
        G = graph(edges(:, 1), edges(:, 2), lengths, nSkel);
    end

    role = repmat("waypoint", nSkel, 1);
    deg = degree(G);
    role(deg == 0) = "island";
    role(deg == 1) = "endpoint";
    role(deg > 2) = "junction";
    role(deg == 2) = "chain";

    skelGraph = struct();
    skelGraph.kind = 'skeleton_graph';
    skelGraph.parent_kind = parentKind;
    skelGraph.node_coords = coords;
    skelGraph.node_labels = fullLabels(fullIdx);
    skelGraph.num_nodes = nSkel;
    skelGraph.edges_local = edges;
    skelGraph.edge_index = edges.';
    skelGraph.edge_length = lengths;
    skelGraph.graph = G;
    skelGraph.boundary_mask = false(nSkel, 1);
    skelGraph.full_node_indices = fullIdx(:);
    skelGraph.full_node_labels = fullLabels(fullIdx);
    skelGraph.grid_linear_index = skelLinear(:);
    skelGraph.boundary_step_distance = distMap(skelLinear);
    skelGraph.node_radius = fullRadius(fullIdx);
    skelGraph.node_thickness = fullThickness(fullIdx);
    skelTable = table( ...
        skelGraph.node_labels(:), ...
        skelGraph.full_node_indices(:), ...
        skelGraph.boundary_step_distance(:), ...
        skelGraph.node_radius(:), ...
        skelGraph.node_thickness(:), ...
        role(:), ...
        'VariableNames', {'Label', 'FullNodeIndex', 'BoundaryStepDistance', 'Radius', 'Thickness', 'SkeletonRole'});
    if istable(fullNodeData) && height(fullNodeData) == size(xyFull, 1)
        attrTable = fullNodeData(fullIdx, :);
        if any(strcmp(attrTable.Properties.VariableNames, 'Label'))
            attrTable(:, 'Label') = [];
        end
        dupNames = intersect(attrTable.Properties.VariableNames, skelTable.Properties.VariableNames, 'stable');
        if ~isempty(dupNames)
            attrTable(:, dupNames) = [];
        end
        skelTable = [skelTable, attrTable];
    end
    skelGraph.node_data_table = skelTable;

    debugData = struct();
    debugData.skelLocalAtGrid = skelLocalAt;
end

function tf = supports_diagonal_link(componentAt, r0, c0, r1, c1, compId)
    tf = true;
    if abs(r1 - r0) ~= 1 || abs(c1 - c0) ~= 1
        return;
    end
    tf = componentAt(r0, c1) == compId || componentAt(r1, c0) == compId;
end

function structGraph = compress_skeleton_graph(skelGraph, detailLevel)
    nSkel = skelGraph.num_nodes;
    if nSkel == 0
        error('extract_structural_graph:EmptySkeletonGraph', 'Skeleton graph has no nodes.');
    end

    if isempty(skelGraph.edges_local)
        structGraph = build_struct_graph_from_segments( ...
            skelGraph, detailLevel, zeros(0, 2), zeros(0, 1), {}, {}, zeros(0, 1), ...
            repmat("island", nSkel, 1), nSkel);
        structGraph.node_coords = skelGraph.node_coords;
        structGraph.node_labels = skelGraph.node_labels(:);
        structGraph.num_nodes = nSkel;
        structGraph.full_node_indices = skelGraph.full_node_indices(:);
        structGraph.full_node_labels = skelGraph.full_node_labels(:);
        structGraph.reduced_node_to_full_indices = skelGraph.full_node_indices(:);
        structGraph.reduced_node_to_skeleton_indices = transpose(1:nSkel);
        structGraph.boundary_mask = false(nSkel, 1);
        if isfield(skelGraph, 'node_radius')
            structGraph.node_radius = skelGraph.node_radius(:);
        end
        if isfield(skelGraph, 'node_thickness')
            structGraph.node_thickness = skelGraph.node_thickness(:);
        end
        structGraph.node_data_table = table( ...
            structGraph.node_labels(:), ...
            structGraph.full_node_indices(:), ...
            repmat("island", nSkel, 1), ...
            'VariableNames', {'Label', 'FullNodeIndex', 'NodeRole'});
        if isfield(skelGraph, 'node_data_table') && istable(skelGraph.node_data_table) ...
                && height(skelGraph.node_data_table) == nSkel
            attrTable = skelGraph.node_data_table;
            if any(strcmp(attrTable.Properties.VariableNames, 'Label'))
                attrTable(:, 'Label') = [];
            end
            dupNames = intersect(attrTable.Properties.VariableNames, structGraph.node_data_table.Properties.VariableNames, 'stable');
            if ~isempty(dupNames)
                attrTable(:, dupNames) = [];
            end
            attrTable.SkeletonNodeIndex = transpose(1:nSkel);
            structGraph.node_data_table = [structGraph.node_data_table, attrTable];
        end
        return;
    end

    G = skelGraph.graph;
    deg = degree(G);
    nEdge = size(skelGraph.edges_local, 1);
    edgeLookup = sparse([skelGraph.edges_local(:, 1); skelGraph.edges_local(:, 2)], ...
        [skelGraph.edges_local(:, 2); skelGraph.edges_local(:, 1)], ...
        [(1:nEdge).'; (1:nEdge).'], ...
        nSkel, nSkel);
    visited = false(size(skelGraph.edges_local, 1), 1);
    compId = conncomp(G);

    reducedNodeForSkel = zeros(nSkel, 1);
    reducedSkelForNode = zeros(0, 1);
    reducedCoords = zeros(0, 2);
    reducedLabels = zeros(0, 1);
    reducedFullIdx = zeros(0, 1);
    reducedRole = strings(0, 1);

    edgePairs = zeros(0, 2);
    edgeLengths = zeros(0, 1);
    edgePolyline = {};
    edgeFullPath = {};
    edgeComp = zeros(0, 1);

    for c = 1:max(compId)
        compNodes = find(compId == c);
        compKey = compNodes(deg(compNodes) ~= 2);

        if isempty(compKey)
            cyclePath = trace_pure_cycle(compNodes, G, edgeLookup);
            add_closed_path(cyclePath, c);
            continue;
        end

        compKey = sort(compKey(:)).';
        for u = compKey
            nbrs = sort(neighbors(G, u)).';
            if isempty(nbrs)
                ensure_reduced_node(u);
                continue;
            end
            for v = nbrs
                eIdx = edgeLookup(u, v);
                if eIdx == 0 || visited(eIdx)
                    continue;
                end
                path = trace_chain(u, v, G, deg, edgeLookup, visited);
                pathEdgeIds = path_to_edge_ids(path, edgeLookup);
                visited(pathEdgeIds) = true;
                if path(1) == path(end)
                    add_closed_path(path, c);
                else
                    add_open_path(path, c);
                end
            end
        end
    end

    freeEdges = find(~visited);
    while ~isempty(freeEdges)
        startEdge = freeEdges(1);
        u = skelGraph.edges_local(startEdge, 1);
        v = skelGraph.edges_local(startEdge, 2);
        cyclePath = trace_chain(u, v, G, deg, edgeLookup, visited);
        pathEdgeIds = path_to_edge_ids(cyclePath, edgeLookup);
        visited(pathEdgeIds) = true;
        if cyclePath(1) == cyclePath(end)
            add_closed_path(cyclePath, compId(u));
        else
            add_open_path(cyclePath, compId(u));
        end
        freeEdges = find(~visited);
    end

    structGraph = build_struct_graph_from_segments( ...
        skelGraph, detailLevel, edgePairs, edgeLengths, edgePolyline, edgeFullPath, edgeComp, ...
        reducedRole, size(reducedCoords, 1));

    structGraph.node_coords = reducedCoords;
    structGraph.node_labels = reducedLabels;
    structGraph.num_nodes = size(reducedCoords, 1);
    structGraph.full_node_indices = reducedFullIdx;
    structGraph.full_node_labels = reducedLabels;
    structGraph.reduced_node_to_full_indices = reducedFullIdx;
    structGraph.reduced_node_to_skeleton_indices = reducedSkelForNode;
    structGraph.boundary_mask = false(structGraph.num_nodes, 1);
    if isfield(skelGraph, 'node_radius')
        structGraph.node_radius = skelGraph.node_radius(reducedSkelForNode);
    end
    if isfield(skelGraph, 'node_thickness')
        structGraph.node_thickness = skelGraph.node_thickness(reducedSkelForNode);
    end
    structTable = table( ...
        structGraph.node_labels(:), ...
        structGraph.full_node_indices(:), ...
        reducedRole(:), ...
        'VariableNames', {'Label', 'FullNodeIndex', 'NodeRole'});
    if isfield(skelGraph, 'node_data_table') && istable(skelGraph.node_data_table) ...
            && height(skelGraph.node_data_table) == skelGraph.num_nodes
        attrTable = skelGraph.node_data_table(reducedSkelForNode, :);
        if any(strcmp(attrTable.Properties.VariableNames, 'Label'))
            attrTable(:, 'Label') = [];
        end
        dupNames = intersect(attrTable.Properties.VariableNames, structTable.Properties.VariableNames, 'stable');
        if ~isempty(dupNames)
            attrTable(:, dupNames) = [];
        end
        attrTable.SkeletonNodeIndex = reducedSkelForNode(:);
        structTable = [structTable, attrTable];
    end
    structGraph.node_data_table = structTable;

    function add_open_path(path, componentId)
        anchorPos = open_path_anchor_positions(path, detailLevel);
        anchorPos = ensure_anchor_positions(anchorPos, numel(path));

        [anchorPos, forced] = ensure_unique_open_path(anchorPos, path);
        if forced
            anchorPos = ensure_anchor_positions(anchorPos, numel(path));
        end

        for segIdx = 1:(numel(anchorPos) - 1)
            segPath = path(anchorPos(segIdx):anchorPos(segIdx + 1));
            add_segment(segPath, componentId);
        end
    end

    function add_closed_path(path, componentId)
        cycleNodes = path(1:end-1);
        anchorPos = closed_path_anchor_positions(cycleNodes, detailLevel);
        anchorPos = unique(anchorPos(:).', 'stable');
        anchorPos = anchorPos(anchorPos >= 1 & anchorPos <= numel(cycleNodes));
        if numel(anchorPos) < 3
            anchorPos = unique([1, round(numel(cycleNodes) / 3), round(2 * numel(cycleNodes) / 3)], 'stable');
            anchorPos = anchorPos(anchorPos >= 1 & anchorPos <= numel(cycleNodes));
        end
        if numel(anchorPos) < 3
            anchorPos = 1:numel(cycleNodes);
        end

        anchorPos = sort(anchorPos);
        wrapNodes = [cycleNodes(:); cycleNodes(anchorPos(1))];
        wrapPos = [anchorPos(:); anchorPos(1) + numel(cycleNodes)];
        for segIdx = 1:(numel(wrapPos) - 1)
            posA = wrapPos(segIdx);
            posB = wrapPos(segIdx + 1);
            segPath = segment_from_cycle(cycleNodes, posA, posB);
            add_segment(segPath, componentId);
        end
    end

    function add_segment(segPath, componentId)
        if numel(segPath) < 2
            return;
        end

        src = ensure_reduced_node(segPath(1));
        dst = ensure_reduced_node(segPath(end));
        if src == dst
            return;
        end

        edgePairs(end + 1, :) = [src, dst]; %#ok<AGROW>
        edgeLengths(end + 1, 1) = path_length(skelGraph.node_coords(segPath, :)); %#ok<AGROW>
        edgePolyline{end + 1, 1} = skelGraph.node_coords(segPath, :); %#ok<AGROW>
        edgeFullPath{end + 1, 1} = skelGraph.full_node_indices(segPath(:)); %#ok<AGROW>
        edgeComp(end + 1, 1) = componentId; %#ok<AGROW>
    end

    function idx = ensure_reduced_node(skelIdx)
        idx = reducedNodeForSkel(skelIdx);
        if idx > 0
            reducedRole(idx) = promote_role(reducedRole(idx), classify_skeleton_node(skelIdx, deg(skelIdx)));
            return;
        end

        idx = size(reducedCoords, 1) + 1;
        reducedNodeForSkel(skelIdx) = idx;
        reducedSkelForNode(idx, 1) = skelIdx; %#ok<AGROW>
        reducedCoords(idx, :) = skelGraph.node_coords(skelIdx, :); %#ok<AGROW>
        reducedLabels(idx, 1) = skelGraph.node_labels(skelIdx); %#ok<AGROW>
        reducedFullIdx(idx, 1) = skelGraph.full_node_indices(skelIdx); %#ok<AGROW>
        reducedRole(idx, 1) = classify_skeleton_node(skelIdx, deg(skelIdx)); %#ok<AGROW>
    end

    function [anchorPos, forced] = ensure_unique_open_path(anchorPos, path)
        forced = false;
        if numel(anchorPos) ~= 2
            return;
        end

        a = reducedNodeForSkel(path(anchorPos(1)));
        b = reducedNodeForSkel(path(anchorPos(2)));
        if a == 0 || b == 0
            return;
        end

        pair = sort([a, b]);
        if isempty(edgePairs)
            return;
        end
        if ~any(all(sort(edgePairs, 2) == pair, 2))
            return;
        end

        internal = 2:(numel(path) - 1);
        if isempty(internal)
            return;
        end

        midPos = internal(round((numel(internal) + 1) / 2));
        anchorPos = [1, midPos, numel(path)];
        forced = true;
    end
end

function path = trace_chain(startNode, nextNode, G, deg, edgeLookup, visited)
    path = [startNode, nextNode];
    prev = startNode;
    curr = nextNode;

    while true
        if curr == startNode && numel(path) > 2
            break;
        end
        if deg(curr) ~= 2
            break;
        end

        nbrs = sort(neighbors(G, curr)).';
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            break;
        end

        next = nbrs(1);
        eIdx = edgeLookup(curr, next);
        if eIdx == 0
            break;
        end
        if visited(eIdx) && next ~= startNode
            break;
        end

        path(end + 1) = next; %#ok<AGROW>
        prev = curr;
        curr = next;
    end
end

function edgeIds = path_to_edge_ids(path, edgeLookup)
    edgeIds = zeros(numel(path) - 1, 1);
    for i = 1:(numel(path) - 1)
        edgeIds(i) = edgeLookup(path(i), path(i + 1));
    end
    edgeIds = edgeIds(edgeIds > 0);
end

function cyclePath = trace_pure_cycle(compNodes, G, edgeLookup)
    compNodes = sort(compNodes(:)).';
    startNode = compNodes(1);
    nbrs = sort(neighbors(G, startNode)).';
    if isempty(nbrs)
        cyclePath = startNode;
        return;
    end

    cyclePath = [startNode, nbrs(1)];
    prev = startNode;
    curr = nbrs(1);

    while true
        nbrs = sort(neighbors(G, curr)).';
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            break;
        end
        next = nbrs(1);
        if edgeLookup(curr, next) == 0
            break;
        end
        cyclePath(end + 1) = next; %#ok<AGROW>
        prev = curr;
        curr = next;
        if curr == startNode
            break;
        end
    end
end

function pos = open_path_anchor_positions(path, detailLevel)
    n = numel(path);
    if n <= 2
        pos = [1, n];
        return;
    end

    nInternal = n - 2;
    nExtra = min(nInternal, round(detailLevel * nInternal));
    if nExtra == 0
        pos = [1, n];
        return;
    end

    internalPos = evenly_spaced_internal_positions(n, nExtra);
    pos = [1, internalPos, n];
end

function pos = closed_path_anchor_positions(cycleNodes, detailLevel)
    n = numel(cycleNodes);
    if n <= 3
        pos = 1:n;
        return;
    end

    nKeep = min(n, max(3, 1 + round(detailLevel * (n - 1))));
    if nKeep >= n
        pos = 1:n;
        return;
    end

    target = linspace(1, n + 1, nKeep + 1);
    pos = unique(min(n, max(1, round(target(1:end-1)))), 'stable');
    if numel(pos) < 3
        pos = unique([1, round(n / 3), round(2 * n / 3)], 'stable');
    end
end

function segPath = segment_from_cycle(cycleNodes, posA, posB)
    n = numel(cycleNodes);
    idx = posA:posB;
    idx = mod(idx - 1, n) + 1;
    segPath = cycleNodes(idx);
end

function pos = evenly_spaced_internal_positions(nNodesInPath, nExtra)
    if nExtra <= 0
        pos = zeros(1, 0);
        return;
    end

    internal = 2:(nNodesInPath - 1);
    target = linspace(2, nNodesInPath - 1, nExtra);
    pos = unique(round(target), 'stable');
    pos = intersect(internal, pos, 'stable');

    k = 1;
    while numel(pos) < nExtra && k <= numel(internal)
        if ~ismember(internal(k), pos)
            pos(end + 1) = internal(k); %#ok<AGROW>
        end
        k = k + 1;
    end

    pos = sort(pos);
end

function pos = ensure_anchor_positions(pos, nPath)
    pos = unique(pos(:).', 'stable');
    pos = pos(pos >= 1 & pos <= nPath);
    pos = unique([1, pos, nPath], 'stable');
end

function out = classify_skeleton_node(skelIdx, deg)
    if deg == 0
        out = "island";
    elseif deg == 1
        out = "endpoint";
    elseif deg > 2
        out = "junction";
    else
        out = "waypoint";
    end
end

function out = promote_role(oldRole, newRole)
    order = ["waypoint", "island", "endpoint", "junction"];
    oldPos = find(order == string(oldRole), 1);
    newPos = find(order == string(newRole), 1);
    if isempty(oldPos)
        oldPos = 1;
    end
    if isempty(newPos)
        newPos = 1;
    end
    out = order(max(oldPos, newPos));
end

function graphOut = build_struct_graph_from_segments(skelGraph, detailLevel, edgePairs, edgeLengths, edgePolyline, edgeFullPath, edgeComp, reducedRole, nNodes)
    graphOut = struct();
    graphOut.kind = 'structural_graph';
    graphOut.parent_kind = 'skeleton_graph';
    graphOut.detail_level = detailLevel;
    graphOut.skeleton_num_nodes = skelGraph.num_nodes;
    graphOut.edge_polyline = edgePolyline;
    graphOut.edge_full_node_paths = edgeFullPath;
    graphOut.edge_component_id = edgeComp;
    graphOut.edge_length = edgeLengths;
    graphOut.edges_local = edgePairs;
    graphOut.edge_index = edgePairs.';
    graphOut.node_role = reducedRole;

    if isempty(edgePairs)
        graphOut.graph = graph([], [], [], nNodes);
    else
        graphOut.graph = graph(edgePairs(:, 1), edgePairs(:, 2), edgeLengths, nNodes);
    end
end

function [graphOut, cleanupDebug] = cleanup_structural_graph(graphIn, skelGraph)
    graphOut = ensure_identity_cleanup_metadata(graphIn);
    cleanupDebug = struct( ...
        'applied', false, ...
        'grid_spacing', NaN, ...
        'min_edge_length', NaN, ...
        'merge_radius', NaN, ...
        'max_iterations', 0, ...
        'iterations_used', 0, ...
        'short_edge_contractions', 0, ...
        'cluster_merges', 0, ...
        'removed_node_coords', zeros(0, size(graphIn.node_coords, 2)), ...
        'removed_original_node_indices', zeros(0, 1), ...
        'before_num_nodes', graphIn.num_nodes, ...
        'after_num_nodes', graphIn.num_nodes, ...
        'before_num_edges', size(graphIn.edges_local, 1), ...
        'after_num_edges', size(graphIn.edges_local, 1));

    if graphIn.num_nodes < 2
        graphOut.cleanup_info = default_cleanup_info(NaN, NaN, NaN, 0, 0, 0);
        return;
    end

    h = estimate_characteristic_spacing(skelGraph.node_coords, skelGraph.edges_local);
    if ~isfinite(h) || h <= 0
        graphOut.cleanup_info = default_cleanup_info(NaN, NaN, NaN, 0, 0, 0);
        return;
    end

    minEdgeLength = 1.5 * h;
    mergeRadius = 1.5 * h;
    maxIterations = 4;

    state = init_cleanup_state(graphIn);
    totalShortEdgeContractions = 0;
    totalClusterMerges = 0;
    removedCoords = zeros(0, size(state.coords, 2));
    removedOriginalIds = zeros(0, 1);
    iterationsUsed = 0;

    for iter = 1:maxIterations
        changed = false;
        iterationsUsed = iter;

        while true
            [edgeIdx, nodeGroup] = find_short_edge_contraction_candidate(state, minEdgeLength);
            if edgeIdx == 0
                break;
            end
            [state, mergeInfo] = merge_node_group(state, nodeGroup);
            if isempty(mergeInfo.removed_indices)
                break;
            end
            removedCoords = [removedCoords; mergeInfo.removed_coords]; %#ok<AGROW>
            removedOriginalIds = [removedOriginalIds; mergeInfo.removed_original_indices(:)]; %#ok<AGROW>
            totalShortEdgeContractions = totalShortEdgeContractions + 1;
            changed = true;
        end

        while true
            nodePair = find_junction_cluster_pair(state, mergeRadius);
            if isempty(nodePair)
                break;
            end
            [state, mergeInfo] = merge_node_group(state, nodePair);
            if isempty(mergeInfo.removed_indices)
                break;
            end
            removedCoords = [removedCoords; mergeInfo.removed_coords]; %#ok<AGROW>
            removedOriginalIds = [removedOriginalIds; mergeInfo.removed_original_indices(:)]; %#ok<AGROW>
            totalClusterMerges = totalClusterMerges + 1;
            changed = true;
        end

        if ~changed
            break;
        end
    end

    if totalShortEdgeContractions == 0 && totalClusterMerges == 0
        cleanupDebug.grid_spacing = h;
        cleanupDebug.min_edge_length = minEdgeLength;
        cleanupDebug.merge_radius = mergeRadius;
        cleanupDebug.max_iterations = maxIterations;
        cleanupDebug.iterations_used = iterationsUsed;
        graphOut.cleanup_info = default_cleanup_info(h, minEdgeLength, mergeRadius, maxIterations, iterationsUsed, 0, 0);
        return;
    end

    graphOut = finalize_cleanup_state(graphIn, state, h, minEdgeLength, mergeRadius, maxIterations, iterationsUsed, ...
        totalShortEdgeContractions, totalClusterMerges);

    cleanupDebug.applied = true;
    cleanupDebug.grid_spacing = h;
    cleanupDebug.min_edge_length = minEdgeLength;
    cleanupDebug.merge_radius = mergeRadius;
    cleanupDebug.max_iterations = maxIterations;
    cleanupDebug.iterations_used = iterationsUsed;
    cleanupDebug.short_edge_contractions = totalShortEdgeContractions;
    cleanupDebug.cluster_merges = totalClusterMerges;
    cleanupDebug.removed_node_coords = removedCoords;
    cleanupDebug.removed_original_node_indices = unique(removedOriginalIds(:), 'stable');
    cleanupDebug.before_num_nodes = graphIn.num_nodes;
    cleanupDebug.after_num_nodes = graphOut.num_nodes;
    cleanupDebug.before_num_edges = size(graphIn.edges_local, 1);
    cleanupDebug.after_num_edges = size(graphOut.edges_local, 1);
end

function graphOut = ensure_identity_cleanup_metadata(graphIn)
    graphOut = graphIn;
    n = graphIn.num_nodes;
    if ~isfield(graphOut, 'merged_full_node_indices') || numel(graphOut.merged_full_node_indices) ~= n
        if isfield(graphOut, 'full_node_indices') && numel(graphOut.full_node_indices) == n
            graphOut.merged_full_node_indices = arrayfun(@(idx) graphOut.full_node_indices(idx), transpose(1:n), 'UniformOutput', false);
        else
            graphOut.merged_full_node_indices = arrayfun(@(idx) graphOut.node_labels(idx), transpose(1:n), 'UniformOutput', false);
        end
    end
    if ~isfield(graphOut, 'merged_skeleton_node_indices') || numel(graphOut.merged_skeleton_node_indices) ~= n
        if isfield(graphOut, 'reduced_node_to_skeleton_indices') && numel(graphOut.reduced_node_to_skeleton_indices) == n
            graphOut.merged_skeleton_node_indices = arrayfun(@(idx) graphOut.reduced_node_to_skeleton_indices(idx), transpose(1:n), 'UniformOutput', false);
        else
            graphOut.merged_skeleton_node_indices = arrayfun(@(idx) idx, transpose(1:n), 'UniformOutput', false);
        end
    end
    if ~isfield(graphOut, 'merged_structural_node_indices') || numel(graphOut.merged_structural_node_indices) ~= n
        graphOut.merged_structural_node_indices = arrayfun(@(idx) idx, transpose(1:n), 'UniformOutput', false);
    end
end

function info = default_cleanup_info(gridSpacing, minEdgeLength, mergeRadius, maxIterations, iterationsUsed, shortEdgeContractions, clusterMerges)
    if nargin < 6
        shortEdgeContractions = 0;
    end
    if nargin < 7
        clusterMerges = 0;
    end
    info = struct( ...
        'grid_spacing', gridSpacing, ...
        'min_edge_length', minEdgeLength, ...
        'merge_radius', mergeRadius, ...
        'max_iterations', maxIterations, ...
        'iterations_used', iterationsUsed, ...
        'short_edge_contractions', shortEdgeContractions, ...
        'cluster_merges', clusterMerges);
end

function state = init_cleanup_state(graphIn)
    n = graphIn.num_nodes;
    state = struct();
    state.coords = graphIn.node_coords;
    state.labels = graphIn.node_labels(:);
    if isfield(graphIn, 'full_node_indices') && numel(graphIn.full_node_indices) == n
        state.full_idx = graphIn.full_node_indices(:);
    else
        state.full_idx = graphIn.node_labels(:);
    end
    if isfield(graphIn, 'reduced_node_to_skeleton_indices') && numel(graphIn.reduced_node_to_skeleton_indices) == n
        state.skel_idx = graphIn.reduced_node_to_skeleton_indices(:);
    else
        state.skel_idx = transpose(1:n);
    end
    if isfield(graphIn, 'node_role') && numel(graphIn.node_role) == n
        state.node_role = string(graphIn.node_role(:));
    else
        state.node_role = repmat("waypoint", n, 1);
    end
    if isfield(graphIn, 'node_radius') && numel(graphIn.node_radius) == n
        state.node_radius = graphIn.node_radius(:);
    else
        state.node_radius = nan(n, 1);
    end
    if isfield(graphIn, 'node_thickness') && numel(graphIn.node_thickness) == n
        state.node_thickness = graphIn.node_thickness(:);
    else
        state.node_thickness = nan(n, 1);
    end
    if isfield(graphIn, 'node_data_table') && istable(graphIn.node_data_table) ...
            && height(graphIn.node_data_table) == n
        state.node_table = graphIn.node_data_table;
    else
        state.node_table = table();
    end

    state.original_struct_idx = arrayfun(@(idx) idx, transpose(1:n), 'UniformOutput', false);
    state.member_full = arrayfun(@(idx) state.full_idx(idx), transpose(1:n), 'UniformOutput', false);
    state.member_skel = arrayfun(@(idx) state.skel_idx(idx), transpose(1:n), 'UniformOutput', false);

    state.edges = graphIn.edges_local;
    nEdge = size(state.edges, 1);
    if isfield(graphIn, 'edge_polyline') && numel(graphIn.edge_polyline) == nEdge
        state.edge_polyline = graphIn.edge_polyline(:);
    else
        state.edge_polyline = arrayfun(@(i) state.coords(state.edges(i, :), :), transpose(1:nEdge), 'UniformOutput', false);
    end
    if isfield(graphIn, 'edge_full_node_paths') && numel(graphIn.edge_full_node_paths) == nEdge
        state.edge_full_path = graphIn.edge_full_node_paths(:);
    else
        state.edge_full_path = arrayfun(@(i) state.full_idx(state.edges(i, :)).', transpose(1:nEdge), 'UniformOutput', false);
    end
    if isfield(graphIn, 'edge_length') && numel(graphIn.edge_length) == nEdge
        state.edge_length = graphIn.edge_length(:);
    else
        state.edge_length = edge_lengths_from_state(state);
    end
end

function h = estimate_characteristic_spacing(coords, edges)
    xy = coords(:, 1:2);
    [xVals, ~] = cluster_axis_values(xy(:, 1));
    [yVals, ~] = cluster_axis_values(xy(:, 2));
    diffs = [diff(sort(xVals(:))); diff(sort(yVals(:)))];
    diffs = diffs(diffs > 0);
    if ~isempty(diffs)
        h = median(diffs);
        return;
    end

    if nargin >= 2 && ~isempty(edges)
        delta = xy(edges(:, 1), :) - xy(edges(:, 2), :);
        edgeLen = sqrt(sum(delta .^ 2, 2));
        edgeLen = edgeLen(edgeLen > 0);
        if ~isempty(edgeLen)
            h = median(edgeLen);
            return;
        end
    end

    h = NaN;
end

function [edgeIdx, nodeGroup] = find_short_edge_contraction_candidate(state, minEdgeLength)
    edgeIdx = 0;
    nodeGroup = zeros(1, 0);
    if isempty(state.edges)
        return;
    end

    G = graph(state.edges(:, 1), state.edges(:, 2), [], size(state.coords, 1));
    deg = degree(G);
    [~, order] = sort(state.edge_length(:), 'ascend');

    for k = 1:numel(order)
        idx = order(k);
        if state.edge_length(idx) >= minEdgeLength
            break;
        end
        u = state.edges(idx, 1);
        v = state.edges(idx, 2);
        if u == v || deg(u) == 1 || deg(v) == 1
            continue;
        end
        if deg(u) ~= 2 || deg(v) ~= 2 || nodes_share_common_neighbor(G, u, v)
            edgeIdx = idx;
            nodeGroup = [u, v];
            return;
        end
    end
end

function tf = nodes_share_common_neighbor(G, u, v)
    nbrU = neighbors(G, u);
    nbrV = neighbors(G, v);
    nbrU(nbrU == v) = [];
    nbrV(nbrV == u) = [];
    tf = ~isempty(intersect(nbrU, nbrV));
end

function nodePair = find_junction_cluster_pair(state, mergeRadius)
    nodePair = zeros(1, 0);
    n = size(state.coords, 1);
    if n < 2
        return;
    end
    if isempty(state.edges)
        return;
    end

    G = graph(state.edges(:, 1), state.edges(:, 2), [], n);
    deg = degree(G);
    compId = conncomp(G);
    candidates = find(deg > 2);
    if numel(candidates) < 2
        return;
    end

    bestDist = inf;
    bestPair = zeros(1, 2);
    xy = state.coords(:, 1:2);
    for i = 1:(numel(candidates) - 1)
        u = candidates(i);
        for j = (i + 1):numel(candidates)
            v = candidates(j);
            if compId(u) ~= compId(v)
                continue;
            end
            d = norm(xy(u, :) - xy(v, :));
            if d <= mergeRadius && d < bestDist
                bestDist = d;
                bestPair = [u, v];
            end
        end
    end

    if all(bestPair > 0)
        nodePair = bestPair;
    end
end

function [state, mergeInfo] = merge_node_group(state, nodeGroup)
    nodeGroup = unique(nodeGroup(:), 'stable');
    nodeGroup = nodeGroup(nodeGroup >= 1 & nodeGroup <= size(state.coords, 1));
    mergeInfo = struct('representative_index', 0, 'removed_indices', zeros(0, 1), ...
        'removed_original_indices', zeros(0, 1), 'removed_coords', zeros(0, size(state.coords, 2)));
    if numel(nodeGroup) < 2
        return;
    end

    G = graph(state.edges(:, 1), state.edges(:, 2), [], size(state.coords, 1));
    deg = degree(G);
    rep = choose_representative_node(state, nodeGroup, deg);
    removeNodes = setdiff(nodeGroup, rep, 'stable');
    if isempty(removeNodes)
        return;
    end

    state.member_full{rep} = combine_unique_numeric_vectors( ...
        state.member_full{rep}, concatenate_numeric_cells(state.member_full(removeNodes)));
    state.member_skel{rep} = combine_unique_numeric_vectors( ...
        state.member_skel{rep}, concatenate_numeric_cells(state.member_skel(removeNodes)));
    state.original_struct_idx{rep} = combine_unique_numeric_vectors( ...
        state.original_struct_idx{rep}, concatenate_numeric_cells(state.original_struct_idx(removeNodes)));

    for e = 1:size(state.edges, 1)
        src = state.edges(e, 1);
        dst = state.edges(e, 2);
        poly = state.edge_polyline{e};

        if ismember(src, removeNodes)
            state.edges(e, 1) = rep;
            if ~isempty(poly)
                poly(1, :) = state.coords(rep, :);
            end
        end
        if ismember(dst, removeNodes)
            state.edges(e, 2) = rep;
            if ~isempty(poly)
                poly(end, :) = state.coords(rep, :);
            end
        end
        state.edge_polyline{e} = poly;
    end

    removeMask = false(size(state.coords, 1), 1);
    removeMask(removeNodes) = true;

    mergeInfo.representative_index = rep;
    mergeInfo.removed_indices = removeNodes(:);
    mergeInfo.removed_original_indices = concatenate_numeric_cells(state.original_struct_idx(removeNodes));
    mergeInfo.removed_coords = state.coords(removeNodes, :);

    state = compact_cleanup_state(state, removeMask);
end

function rep = choose_representative_node(state, nodeGroup, deg)
    degVals = deg(nodeGroup);
    bestDegree = max(degVals);
    candidates = nodeGroup(degVals == bestDegree);

    thicknessVals = node_thickness_proxy(state);
    if ~isempty(thicknessVals)
        thick = thicknessVals(candidates);
        thick(~isfinite(thick)) = -inf;
        bestThickness = max(thick);
        candidates = candidates(thick == bestThickness);
    end

    if numel(candidates) == 1
        rep = candidates(1);
        return;
    end

    xy = state.coords(candidates, 1:2);
    totalDist = zeros(numel(candidates), 1);
    for i = 1:numel(candidates)
        delta = xy - xy(i, :);
        totalDist(i) = sum(sqrt(sum(delta .^ 2, 2)));
    end
    [~, bestIdx] = min(totalDist);
    rep = candidates(bestIdx);
end

function vals = node_thickness_proxy(state)
    vals = [];
    if isfield(state, 'node_thickness') && numel(state.node_thickness) == size(state.coords, 1)
        vals = double(state.node_thickness(:));
        return;
    end
    if ~istable(state.node_table) || isempty(state.node_table)
        return;
    end

    preferred = {'Thickness', 'LocalThickness', 'BoundaryStepDistance', ...
        'Attr_Thickness', 'Attr_LocalThickness', 'Attr_BoundaryStepDistance'};
    varNames = state.node_table.Properties.VariableNames;
    for i = 1:numel(preferred)
        if any(strcmp(varNames, preferred{i}))
            col = state.node_table.(preferred{i});
            if isnumeric(col) && numel(col) == size(state.coords, 1)
                vals = double(col(:));
                return;
            end
        end
    end
end

function out = concatenate_numeric_cells(cellVals)
    out = zeros(0, 1);
    for i = 1:numel(cellVals)
        vals = cellVals{i};
        if isempty(vals)
            continue;
        end
        out = [out; vals(:)]; %#ok<AGROW>
    end
    out = unique(out, 'stable');
end

function out = combine_unique_numeric_vectors(a, b)
    out = unique([a(:); b(:)], 'stable');
end

function state = compact_cleanup_state(state, removeMask)
    keepMask = ~removeMask(:);
    oldToNew = zeros(numel(keepMask), 1);
    oldToNew(keepMask) = 1:nnz(keepMask);

    state.coords = state.coords(keepMask, :);
    state.labels = state.labels(keepMask);
    state.full_idx = state.full_idx(keepMask);
    state.skel_idx = state.skel_idx(keepMask);
    state.node_role = state.node_role(keepMask);
    state.node_radius = state.node_radius(keepMask);
    state.node_thickness = state.node_thickness(keepMask);
    state.member_full = state.member_full(keepMask);
    state.member_skel = state.member_skel(keepMask);
    state.original_struct_idx = state.original_struct_idx(keepMask);
    if istable(state.node_table) && height(state.node_table) == numel(keepMask)
        state.node_table = state.node_table(keepMask, :);
    end

    if isempty(state.edges)
        return;
    end

    state.edges = oldToNew(state.edges);
    valid = all(state.edges > 0, 2);
    state.edges = state.edges(valid, :);
    state.edge_polyline = state.edge_polyline(valid);
    state.edge_full_path = state.edge_full_path(valid);
    state.edge_length = state.edge_length(valid);

    nonSelf = state.edges(:, 1) ~= state.edges(:, 2);
    state.edges = state.edges(nonSelf, :);
    state.edge_polyline = state.edge_polyline(nonSelf);
    state.edge_full_path = state.edge_full_path(nonSelf);
    state.edge_length = state.edge_length(nonSelf);

    state = deduplicate_cleanup_edges(state);
end

function state = deduplicate_cleanup_edges(state)
    if isempty(state.edges)
        return;
    end

    for i = 1:size(state.edges, 1)
        state.edge_length(i, 1) = resolve_edge_length(state, i);
    end

    sortedPairs = sort(state.edges, 2);
    [~, ~, groupId] = unique(sortedPairs, 'rows', 'stable');
    keepIdx = zeros(max(groupId), 1);
    for g = 1:max(groupId)
        members = find(groupId == g);
        if numel(members) == 1
            keepIdx(g) = members;
            continue;
        end
        [~, bestRel] = max(state.edge_length(members));
        keepIdx(g) = members(bestRel(1));
    end
    keepIdx = keepIdx(keepIdx > 0);
    keepIdx = sort(keepIdx);

    state.edges = state.edges(keepIdx, :);
    state.edge_polyline = state.edge_polyline(keepIdx);
    state.edge_full_path = state.edge_full_path(keepIdx);
    state.edge_length = state.edge_length(keepIdx);
end

function len = resolve_edge_length(state, edgeIdx)
    poly = state.edge_polyline{edgeIdx};
    if ~isempty(poly) && size(poly, 1) >= 2
        len = path_length(poly);
        return;
    end
    edge = state.edges(edgeIdx, :);
    len = norm(state.coords(edge(1), :) - state.coords(edge(2), :));
end

function lengths = edge_lengths_from_state(state)
    lengths = zeros(size(state.edges, 1), 1);
    for i = 1:size(state.edges, 1)
        lengths(i, 1) = resolve_edge_length(state, i);
    end
end

function graphOut = finalize_cleanup_state(graphIn, state, h, minEdgeLength, mergeRadius, maxIterations, iterationsUsed, shortEdgeContractions, clusterMerges)
    graphOut = graphIn;
    n = size(state.coords, 1);
    graphOut.node_coords = state.coords;
    graphOut.node_labels = state.labels;
    graphOut.num_nodes = n;
    graphOut.full_node_indices = state.full_idx;
    graphOut.full_node_labels = state.labels;
    graphOut.reduced_node_to_full_indices = state.full_idx;
    graphOut.reduced_node_to_skeleton_indices = state.skel_idx;
    graphOut.boundary_mask = false(n, 1);
    graphOut.node_radius = state.node_radius;
    graphOut.node_thickness = state.node_thickness;
    graphOut.merged_full_node_indices = state.member_full;
    graphOut.merged_skeleton_node_indices = state.member_skel;
    graphOut.merged_structural_node_indices = state.original_struct_idx;

    graphOut.edges_local = state.edges;
    graphOut.edge_index = state.edges.';
    graphOut.edge_polyline = state.edge_polyline;
    graphOut.edge_full_node_paths = state.edge_full_path;
    graphOut.edge_length = edge_lengths_from_state(state);

    if isempty(state.edges)
        G = graph([], [], [], n);
        edgeComp = zeros(0, 1);
        deg = zeros(n, 1);
    else
        G = graph(state.edges(:, 1), state.edges(:, 2), graphOut.edge_length, n);
        compId = conncomp(G);
        edgeComp = compId(state.edges(:, 1)).';
        edgeComp = edgeComp(:);
        deg = degree(G);
    end

    graphOut.graph = G;
    graphOut.edge_component_id = edgeComp;
    roleCell = arrayfun(@(d) classify_skeleton_node(1, d), deg, 'UniformOutput', false);
    graphOut.node_role = reshape(string(roleCell), [], 1);

    nodeTable = state.node_table;
    if ~istable(nodeTable) || height(nodeTable) ~= n
        nodeTable = table();
    else
        if any(strcmp(nodeTable.Properties.VariableNames, 'Label'))
            nodeTable.Label = state.labels;
        end
    end

    nodeTable = upsert_table_column(nodeTable, 'FullNodeIndex', state.full_idx);
    nodeTable = upsert_table_column(nodeTable, 'SkeletonNodeIndex', state.skel_idx);
    nodeTable = upsert_table_column(nodeTable, 'Radius', state.node_radius);
    nodeTable = upsert_table_column(nodeTable, 'Thickness', state.node_thickness);
    nodeTable = upsert_table_column(nodeTable, 'NodeRole', graphOut.node_role);
    nodeTable = upsert_table_column(nodeTable, 'MergedFullNodeIndices', reshape(string(cellfun(@join_numeric_path_local, state.member_full, 'UniformOutput', false)), [], 1));
    nodeTable = upsert_table_column(nodeTable, 'MergedSkeletonNodeIndices', reshape(string(cellfun(@join_numeric_path_local, state.member_skel, 'UniformOutput', false)), [], 1));
    nodeTable = upsert_table_column(nodeTable, 'MergedStructuralNodeIndices', reshape(string(cellfun(@join_numeric_path_local, state.original_struct_idx, 'UniformOutput', false)), [], 1));
    graphOut.node_data_table = nodeTable;

    graphOut.cleanup_info = default_cleanup_info(h, minEdgeLength, mergeRadius, maxIterations, iterationsUsed, shortEdgeContractions, clusterMerges);
end

function graphOut = attach_structural_edge_thickness(graphIn, fullGraph)
    graphOut = graphIn;
    nEdge = size(graphIn.edges_local, 1);
    if nEdge == 0 || ~isfield(graphIn, 'edge_full_node_paths') || numel(graphIn.edge_full_node_paths) ~= nEdge ...
            || ~isfield(fullGraph, 'node_thickness') || isempty(fullGraph.node_thickness)
        graphOut.edge_thickness_min = zeros(0, 1);
        graphOut.edge_thickness_mean = zeros(0, 1);
        graphOut.edge_thickness_max = zeros(0, 1);
        graphOut.edge_thickness_bottleneck = zeros(0, 1);
        graphOut.edge_thickness_delta = zeros(0, 1);
        return;
    end

    fullThickness = fullGraph.node_thickness(:);
    edgeMin = nan(nEdge, 1);
    edgeMean = nan(nEdge, 1);
    edgeMax = nan(nEdge, 1);
    edgeBottleneck = nan(nEdge, 1);
    edgeDelta = nan(nEdge, 1);
    for i = 1:nEdge
        pathIdx = graphIn.edge_full_node_paths{i};
        if isempty(pathIdx)
            continue;
        end
        pathIdx = pathIdx(:);
        valid = pathIdx >= 1 & pathIdx <= numel(fullThickness);
        if ~any(valid)
            continue;
        end
        values = fullThickness(pathIdx(valid));
        values = values(isfinite(values));
        if isempty(values)
            continue;
        end
        edgeMin(i) = min(values);
        edgeMean(i) = mean(values);
        edgeMax(i) = max(values);
        edgeBottleneck(i) = edgeMin(i);
        edgeDelta(i) = edgeMax(i) - edgeMin(i);
    end

    graphOut.edge_thickness_min = edgeMin;
    graphOut.edge_thickness_mean = edgeMean;
    graphOut.edge_thickness_max = edgeMax;
    graphOut.edge_thickness_bottleneck = edgeBottleneck;
    graphOut.edge_thickness_delta = edgeDelta;
end

function tbl = upsert_table_column(tbl, varName, values)
    if ~istable(tbl)
        tbl = table();
    end
    if isempty(tbl)
        tbl = table(values, 'VariableNames', {varName});
        return;
    end
    if any(strcmp(tbl.Properties.VariableNames, varName))
        tbl.(varName) = values;
    else
        tbl.(varName) = values;
    end
end

function out = join_numeric_path_local(pathVec)
    if isempty(pathVec)
        out = "";
        return;
    end
    out = strjoin(string(pathVec(:).'), ' ');
end

function out = get_default(s, fieldName, fallback)
    if isstruct(s) && isfield(s, fieldName)
        out = s.(fieldName);
    else
        out = fallback;
    end
end

function len = path_length(coords)
    if size(coords, 1) < 2
        len = 0;
        return;
    end
    delta = diff(coords, 1, 1);
    len = sum(sqrt(sum(delta .^ 2, 2)));
end
