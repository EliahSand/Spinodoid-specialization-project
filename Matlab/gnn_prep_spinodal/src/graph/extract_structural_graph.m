function [structGraph, skeletonGraph, debugData, annotatedFullGraph] = extract_structural_graph(fullGraph, varargin)
%EXTRACT_STRUCTURAL_GRAPH Build a simplified skeleton graph from a dense graph.
%
% [structGraph, skeletonGraph, debugData, annotatedFullGraph] = EXTRACT_STRUCTURAL_GRAPH(fullGraph, ...)
%
% Name-Value options:
%   'DetailLevel'           : [0,1], 0 = coarsest, 1 = keep full skeleton chains
%   'MinIslandNodes'        : minimum skeleton-node count kept for an isolated component
%   'SpurPruneRadiusFactor' : prune terminal branches shorter than this multiple of junction radius

    p = inputParser;
    p.addParameter('DetailLevel', 0.20, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('MinIslandNodes', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('SpurPruneRadiusFactor', 1.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    p.parse(varargin{:});
    opts = p.Results;

    for f = {'node_coords', 'node_labels', 'edges_local', 'boundary_mask'}
        if ~isfield(fullGraph, f{1})
            error('extract_structural_graph:MissingField', ...
                'fullGraph is missing required field "%s".', f{1});
        end
    end

    xy = fullGraph.node_coords(:, 1:2);
    nFull = size(xy, 1);

    [xVals, xBin] = cluster_axis_values(xy(:, 1));
    [yVals, yBin] = cluster_axis_values(xy(:, 2));
    nRows = numel(yVals);
    nCols = numel(xVals);

    nodeAt = zeros(nRows, nCols);
    for i = 1:nFull
        if nodeAt(yBin(i), xBin(i)) == 0
            nodeAt(yBin(i), xBin(i)) = i;
        end
    end
    occMask = nodeAt > 0;
    componentAt = periodic_components(occMask);

    % Periodic boundary: occupied iff some 4-neighbour (with wrap) is void.
    boundaryMask = occMask & ~( ...
        circshift(occMask,  1, 1) & circshift(occMask, -1, 1) & ...
        circshift(occMask,  1, 2) & circshift(occMask, -1, 2));

    boundaryStepDistance = multi_source_grid_distance(occMask, boundaryMask);
    [radiusMap, thicknessMap] = compute_physical_radius_map(occMask, xVals, yVals);
    annotatedFullGraph = attach_full_graph_thickness(fullGraph, xBin, yBin, radiusMap, thicknessMap);

    % Skeletonize the real occupied material only. Closing the mask before
    % skeletonization can erase local medial branches and create false
    % topological bridges through one-cell void gaps.
    fullSkel = skeletonize(occMask);

    skelMask = false(size(occMask));
    nComp = max(componentAt(:));
    for compId = 1:nComp
        compMask = componentAt == compId;
        if nnz(compMask) < opts.MinIslandNodes
            continue;
        end
        compSkel = fullSkel & compMask;
        % Ensure every 8-connected sub-piece of this periodic component (a
        % single periodic component can split into several pieces when
        % unwrapped) has at least one skeleton representative.
        sub = bwlabel(compMask, 8);
        for subId = 1:max(sub(:))
            pieceMask = sub == subId;
            if nnz(compSkel & pieceMask) >= opts.MinIslandNodes
                continue;
            end
            pieceIdx = find(pieceMask);
            [~, pickRel] = max(boundaryStepDistance(pieceIdx));
            compSkel(pieceIdx(pickRel)) = true;
        end
        skelMask = skelMask | compSkel;
    end

    skeletonGraph = build_skeleton_graph(skelMask, componentAt, nodeAt, xy, ...
        annotatedFullGraph.node_labels(:), annotatedFullGraph.node_radius(:), annotatedFullGraph.node_thickness(:));
    skeletonGraph = prune_radius_covered_spurs(skeletonGraph, opts.SpurPruneRadiusFactor);

    structGraph = compress_skeleton(skeletonGraph, opts.DetailLevel);
    structGraphBeforeCleanup = structGraph;
    [structGraph, cleanupInfo] = cleanup_structural(structGraph);
    structGraph.cleanup_info = cleanupInfo;
    structGraph = attach_edge_thickness_min(structGraph, annotatedFullGraph.node_thickness(:));

    debugData = struct();
    debugData.structural_before_cleanup = structGraphBeforeCleanup;
end

% =========================================================================
% Grid / mask helpers
% =========================================================================

function [valsOut, binIdx] = cluster_axis_values(vals)
    vals = double(vals(:));
    n = numel(vals);
    if n == 0
        valsOut = zeros(0, 1);
        binIdx  = zeros(0, 1);
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

    occupiedVals    = accumarray(groupSorted, sortedVals, [], @mean);
    occupiedBinIdx  = zeros(n, 1);
    occupiedBinIdx(order) = groupSorted;

    [valsOut, occupiedToFull] = expand_regular_axis_values(occupiedVals, tol);
    binIdx = occupiedToFull(occupiedBinIdx);
end

function [valsOut, occupiedToFull] = expand_regular_axis_values(occupiedVals, clusterTol)
    occupiedVals = double(occupiedVals(:));
    n = numel(occupiedVals);
    occupiedToFull = transpose(1:n);
    valsOut = occupiedVals;
    if n < 2
        return;
    end

    diffs = diff(occupiedVals);
    posDiffs = diffs(diffs > 0);
    if isempty(posDiffs)
        return;
    end
    baseStep = min(posDiffs);
    if ~isfinite(baseStep) || baseStep <= 0
        return;
    end

    relPos = (occupiedVals - occupiedVals(1)) ./ baseStep;
    fullIdx = round(relPos) + 1;
    axisTol = max([0.20 * baseStep, 10 * clusterTol, eps(max(abs(occupiedVals))) * 16]);
    reconstructed = occupiedVals(1) + (fullIdx - 1) .* baseStep;
    if any(fullIdx < 1) || any(abs(reconstructed - occupiedVals) > axisTol)
        return;
    end

    nFull = fullIdx(end);
    if nFull <= n || nFull > max(4 * n, n + 256)
        return;
    end

    valsOut = occupiedVals(1) + transpose(0:(nFull - 1)) .* baseStep;
    occupiedToFull = fullIdx(:);
end

function labels = periodic_components(mask)
    [nRows, nCols] = size(mask);
    labels = zeros(nRows, nCols);
    if ~any(mask(:))
        return;
    end
    tiled = repmat(logical(mask), 3, 3);
    tiledLbl = bwlabel(tiled, 8);
    centre = tiledLbl((nRows + 1):(2 * nRows), (nCols + 1):(2 * nCols));
    centre(~mask) = 0;
    [u, ~, ic] = unique(centre(:));
    remap = zeros(numel(u), 1);
    keep = u > 0;
    remap(keep) = 1:nnz(keep);
    labels = reshape(remap(ic), nRows, nCols);
end

function distMap = multi_source_grid_distance(mask, seedMask)
    [nRows, nCols] = size(mask);
    distMap = inf(nRows, nCols);
    seedMask = logical(seedMask) & logical(mask);

    [seedRows, seedCols] = find(seedMask);
    nSeeds = numel(seedRows);
    queue = zeros(nnz(mask), 2);
    queue(1:nSeeds, :) = [seedRows, seedCols];
    if nSeeds > 0
        distMap(sub2ind([nRows, nCols], seedRows, seedCols)) = 0;
    end
    qHead = 1;
    qTail = nSeeds;
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
            if base + 1 < distMap(rr, cc)
                distMap(rr, cc) = base + 1;
                qTail = qTail + 1;
                queue(qTail, :) = [rr, cc];
            end
        end
    end
end

function [radiusMap, thicknessMap] = compute_physical_radius_map(mask, xVals, yVals)
    [nRows, nCols] = size(mask);
    radiusMap    = zeros(nRows, nCols);
    thicknessMap = zeros(nRows, nCols);

    dx = diff(sort(xVals(:)));
    dy = diff(sort(yVals(:)));
    allDiffs = [dx(dx > 0); dy(dy > 0)];
    allDiffs = allDiffs(isfinite(allDiffs));
    if isempty(allDiffs)
        return;
    end
    h = median(allDiffs);
    if ~any(mask(:)) || ~isfinite(h) || h <= 0
        return;
    end

    % Non-periodic Euclidean distance: material genuinely terminates at the
    % domain edges (clamped/loaded/free in the FEA), so radii must reflect
    % that. Tiling would falsely inflate edge-node radii.
    distPix = bwdist(~mask, 'euclidean');

    radiusMap(mask)    = double(distPix(mask)) * h;
    thicknessMap(mask) = 2 * radiusMap(mask);
end

function graphOut = attach_full_graph_thickness(graphIn, xBin, yBin, radiusMap, thicknessMap)
    graphOut = graphIn;
    n = numel(graphIn.node_labels);
    radius    = zeros(n, 1);
    thickness = zeros(n, 1);
    for i = 1:n
        radius(i)    = radiusMap(yBin(i), xBin(i));
        thickness(i) = thicknessMap(yBin(i), xBin(i));
    end
    graphOut.node_radius    = radius;
    graphOut.node_thickness = thickness;

    if ~isfield(graphOut, 'node_data_table') || ~istable(graphOut.node_data_table) ...
            || height(graphOut.node_data_table) ~= n
        graphOut.node_data_table = table(graphOut.node_labels(:), 'VariableNames', {'Label'});
    end
    graphOut.node_data_table.Radius    = radius;
    graphOut.node_data_table.Thickness = thickness;
end

function skel = skeletonize(mask)
    mask = logical(mask);
    if ~any(mask(:))
        skel = false(size(mask));
        return;
    end
    try
        skel = bwmorph(mask, 'skel', Inf);
    catch
        skel = bwskel(mask);
    end
    skel = logical(skel);
end

% =========================================================================
% Skeleton graph
% =========================================================================

function skelGraph = build_skeleton_graph(skelMask, componentAt, nodeAt, xyFull, fullLabels, fullRadius, fullThickness)
    skelLinear = find(skelMask);
    nSkel = numel(skelLinear);
    if nSkel == 0
        error('extract_structural_graph:EmptySkeleton', 'Skeleton extraction produced no nodes.');
    end

    fullIdx = nodeAt(skelLinear);
    coords  = xyFull(fullIdx, :);

    skelLocalAt = zeros(size(skelMask));
    skelLocalAt(skelLinear) = 1:nSkel;

    % Half-circle 8-connected neighbours so each edge is emitted once.
    dirs = [0 1; 1 -1; 1 0; 1 1];
    [nRows, nCols] = size(skelMask);
    maxEdges = nSkel * size(dirs, 1);
    edges = zeros(maxEdges, 2);
    lengths = zeros(maxEdges, 1);
    nE = 0;
    for i = 1:nSkel
        [r, c] = ind2sub([nRows, nCols], skelLinear(i));
        thisComp = componentAt(r, c);
        if thisComp == 0
            continue;
        end
        for d = 1:size(dirs, 1)
            rr = r + dirs(d, 1);
            cc = c + dirs(d, 2);
            if rr < 1 || rr > nRows || cc < 1 || cc > nCols
                continue;
            end
            if componentAt(rr, cc) ~= thisComp
                continue;
            end
            if abs(dirs(d, 1)) == 1 && abs(dirs(d, 2)) == 1
                % Avoid tiny triangles when an orthogonal skeleton pixel
                % already connects this corner, and avoid corner-touch
                % shortcuts across distinct components.
                onOrthSkel = (componentAt(r, cc) == thisComp && skelLocalAt(r, cc) > 0) || ...
                             (componentAt(rr, c) == thisComp && skelLocalAt(rr, c) > 0);
                if onOrthSkel
                    continue;
                end
                if componentAt(r, cc) ~= thisComp && componentAt(rr, c) ~= thisComp
                    continue;
                end
            end
            j = skelLocalAt(rr, cc);
            if j == 0
                continue;
            end
            nE = nE + 1;
            edges(nE, :) = [i, j];
            lengths(nE) = norm(coords(i, :) - coords(j, :));
        end
    end
    edges = edges(1:nE, :);
    lengths = lengths(1:nE);

    skelGraph = struct();
    skelGraph.node_coords    = coords;
    skelGraph.node_labels    = fullLabels(fullIdx);
    skelGraph.num_nodes      = nSkel;
    skelGraph.edges_local    = edges;
    skelGraph.edge_length    = lengths;
    skelGraph.full_node_indices = fullIdx(:);
    skelGraph.node_radius    = fullRadius(fullIdx);
    skelGraph.node_thickness = fullThickness(fullIdx);
    skelGraph.node_data_table = table( ...
        skelGraph.node_labels(:), ...
        skelGraph.full_node_indices, ...
        skelGraph.node_radius, ...
        skelGraph.node_thickness, ...
        'VariableNames', {'Label', 'FullNodeIndex', 'Radius', 'Thickness'});
end

function skelGraph = prune_radius_covered_spurs(skelGraph, radiusFactor)
    if radiusFactor <= 0 || skelGraph.num_nodes < 2 || isempty(skelGraph.edges_local)
        return;
    end

    n = skelGraph.num_nodes;
    G = graph(skelGraph.edges_local(:, 1), skelGraph.edges_local(:, 2), [], n);
    deg = degree(G);
    endpoints = find(deg == 1);
    removeMask = false(n, 1);

    for i = 1:numel(endpoints)
        [path, root] = trace_terminal_branch(G, deg, endpoints(i));
        if isempty(path) || root < 1 || deg(root) <= 2
            continue;
        end
        rootRadius = skelGraph.node_radius(root);
        if ~isfinite(rootRadius) || rootRadius <= 0
            continue;
        end
        if path_length(skelGraph.node_coords(path, 1:2)) <= radiusFactor * rootRadius
            removeMask(path(1:end - 1)) = true;
        end
    end

    if any(removeMask)
        skelGraph = remove_skeleton_nodes(skelGraph, removeMask);
    end
end

function [path, root] = trace_terminal_branch(G, deg, endpoint)
    path = endpoint;
    prev = 0;
    curr = endpoint;
    while true
        nbrs = neighbors(G, curr);
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            root = curr;
            return;
        end
        next = nbrs(1);
        path(end + 1) = next; %#ok<AGROW>
        prev = curr;
        curr = next;
        root = curr;
        if deg(curr) ~= 2
            return;
        end
    end
end

function skelGraph = remove_skeleton_nodes(skelGraph, removeMask)
    keepMask = ~removeMask(:);
    oldToNew = zeros(numel(keepMask), 1);
    oldToNew(keepMask) = 1:nnz(keepMask);

    skelGraph.node_coords       = skelGraph.node_coords(keepMask, :);
    skelGraph.node_labels       = skelGraph.node_labels(keepMask);
    skelGraph.num_nodes         = nnz(keepMask);
    skelGraph.full_node_indices = skelGraph.full_node_indices(keepMask);
    skelGraph.node_radius       = skelGraph.node_radius(keepMask);
    skelGraph.node_thickness    = skelGraph.node_thickness(keepMask);
    skelGraph.node_data_table   = skelGraph.node_data_table(keepMask, :);

    if isempty(skelGraph.edges_local)
        return;
    end
    edges = skelGraph.edges_local;
    edgeKeep = ~removeMask(edges(:, 1)) & ~removeMask(edges(:, 2));
    edges = oldToNew(edges(edgeKeep, :));
    edges = edges(edges(:, 1) ~= edges(:, 2), :);
    skelGraph.edges_local = edges;
    if isempty(edges)
        skelGraph.edge_length = zeros(0, 1);
    else
        delta = skelGraph.node_coords(edges(:, 1), 1:2) - skelGraph.node_coords(edges(:, 2), 1:2);
        skelGraph.edge_length = sqrt(sum(delta .^ 2, 2));
    end
end

% =========================================================================
% Compress chains -> structural graph
% =========================================================================

function structGraph = compress_skeleton(skelGraph, detailLevel)
    nSkel = skelGraph.num_nodes;
    if nSkel == 0
        error('extract_structural_graph:EmptySkeletonGraph', 'Skeleton graph has no nodes.');
    end

    % Shared mutable state (nested fns mutate these via shared workspace).
    reducedForSkel = zeros(nSkel, 1);
    skelForReduced = zeros(0, 1);
    edgePairs    = zeros(0, 2);
    edgeLengths  = zeros(0, 1);
    edgePolyline = {};
    edgeFullPath = {};

    if isempty(skelGraph.edges_local)
        % All-islands case: every skeleton node becomes a structural node.
        skelForReduced = transpose(1:nSkel);
        structGraph = finalize_structural(skelGraph, skelForReduced, ...
            zeros(0, 2), zeros(0, 1), {}, {});
        return;
    end

    G = graph(skelGraph.edges_local(:, 1), skelGraph.edges_local(:, 2), [], nSkel);
    deg = degree(G);
    nE = size(skelGraph.edges_local, 1);
    edgeLookup = sparse([skelGraph.edges_local(:, 1); skelGraph.edges_local(:, 2)], ...
                        [skelGraph.edges_local(:, 2); skelGraph.edges_local(:, 1)], ...
                        [(1:nE).'; (1:nE).'], nSkel, nSkel);
    visited = false(nE, 1);
    compId = conncomp(G);

    for c = 1:max(compId)
        compNodes = find(compId == c);
        anchors = compNodes(deg(compNodes) ~= 2);

        if isempty(anchors)
            cyclePath = trace_pure_cycle(compNodes, G);
            emit_cycle(cyclePath);
            continue;
        end

        anchors = sort(anchors(:)).';
        for u = anchors
            nbrs = sort(neighbors(G, u)).';
            for v = nbrs
                eIdx = edgeLookup(u, v);
                if eIdx == 0 || visited(eIdx)
                    continue;
                end
                path = trace_chain(u, v, G, deg, edgeLookup, visited);
                pathEdges = path_to_edge_ids(path, edgeLookup);
                visited(pathEdges) = true;
                if path(1) == path(end)
                    emit_cycle(path);
                else
                    emit_open(path);
                end
            end
        end
    end

    % Any remaining unvisited edges form free chains (rare/degenerate).
    freeEdges = find(~visited);
    while ~isempty(freeEdges)
        startEdge = freeEdges(1);
        u = skelGraph.edges_local(startEdge, 1);
        v = skelGraph.edges_local(startEdge, 2);
        path = trace_chain(u, v, G, deg, edgeLookup, visited);
        pathEdges = path_to_edge_ids(path, edgeLookup);
        visited(pathEdges) = true;
        if path(1) == path(end)
            emit_cycle(path);
        else
            emit_open(path);
        end
        freeEdges = find(~visited);
    end

    structGraph = finalize_structural(skelGraph, skelForReduced, ...
        edgePairs, edgeLengths, edgePolyline, edgeFullPath);

    % --- nested functions: share workspace with compress_skeleton -------------

    function rIdx = ensure_reduced(sIdx)
        rIdx = reducedForSkel(sIdx);
        if rIdx == 0
            rIdx = numel(skelForReduced) + 1;
            reducedForSkel(sIdx) = rIdx;
            skelForReduced(rIdx, 1) = sIdx;
        end
    end

    function emit_segment(segPath)
        if numel(segPath) < 2
            return;
        end
        src = ensure_reduced(segPath(1));
        dst = ensure_reduced(segPath(end));
        if src == dst
            return;
        end
        edgePairs(end + 1, :)    = [src, dst];
        edgeLengths(end + 1, 1)  = path_length(skelGraph.node_coords(segPath, :));
        edgePolyline{end + 1, 1} = skelGraph.node_coords(segPath, :);
        edgeFullPath{end + 1, 1} = skelGraph.full_node_indices(segPath(:));
    end

    function emit_open(path)
        n = numel(path);
        nInternal = max(n - 2, 0);
        nExtra = min(nInternal, round(detailLevel * nInternal));
        if nExtra == 0
            anchorsLocal = [1, n];
        else
            target = linspace(2, n - 1, nExtra);
            internal = unique(round(target), 'stable');
            internal = internal(internal >= 2 & internal <= n - 1);
            anchorsLocal = [1, internal, n];
        end
        anchorsLocal = unique(anchorsLocal, 'stable');

        % Avoid parallel-edge duplication between two existing reduced nodes.
        if numel(anchorsLocal) == 2 && ~isempty(edgePairs)
            a = reducedForSkel(path(anchorsLocal(1)));
            b = reducedForSkel(path(anchorsLocal(2)));
            if a > 0 && b > 0 && any(all(sort(edgePairs, 2) == sort([a, b]), 2))
                internal = 2:(n - 1);
                if ~isempty(internal)
                    midPos = internal(round((numel(internal) + 1) / 2));
                    anchorsLocal = [1, midPos, n];
                end
            end
        end

        for s = 1:(numel(anchorsLocal) - 1)
            emit_segment(path(anchorsLocal(s):anchorsLocal(s + 1)));
        end
    end

    function emit_cycle(path)
        if numel(path) > 1 && path(1) == path(end)
            cycleNodes = path(1:end - 1);
        else
            cycleNodes = path;
        end
        n = numel(cycleNodes);
        if n == 0
            return;
        end
        if n <= 3
            anchorsLocal = 1:n;
        else
            nKeep = min(n, max(3, 1 + round(detailLevel * (n - 1))));
            if nKeep >= n
                anchorsLocal = 1:n;
            else
                target = linspace(1, n + 1, nKeep + 1);
                anchorsLocal = unique(min(n, max(1, round(target(1:end - 1)))), 'stable');
                if numel(anchorsLocal) < 3
                    anchorsLocal = unique([1, round(n / 3), round(2 * n / 3)], 'stable');
                end
            end
        end
        anchorsLocal = sort(anchorsLocal);

        wrapPos = [anchorsLocal(:); anchorsLocal(1) + n];
        for s = 1:(numel(wrapPos) - 1)
            idx = mod((wrapPos(s):wrapPos(s + 1)) - 1, n) + 1;
            emit_segment(cycleNodes(idx));
        end
    end
end

function structGraph = finalize_structural(skelGraph, skelForReduced, edgePairs, edgeLengths, edgePolyline, edgeFullPath)
    n = numel(skelForReduced);
    structGraph = struct();
    structGraph.node_coords    = skelGraph.node_coords(skelForReduced, :);
    structGraph.node_labels    = skelGraph.node_labels(skelForReduced);
    structGraph.num_nodes      = n;
    structGraph.edges_local    = edgePairs;
    structGraph.edge_index     = edgePairs.';
    structGraph.edge_length    = edgeLengths;
    structGraph.edge_polyline  = edgePolyline;
    structGraph.edge_full_node_paths = edgeFullPath;
    structGraph.full_node_indices = skelGraph.full_node_indices(skelForReduced);
    structGraph.reduced_node_to_full_indices = structGraph.full_node_indices;
    structGraph.node_radius    = skelGraph.node_radius(skelForReduced);
    structGraph.node_thickness = skelGraph.node_thickness(skelForReduced);
    structGraph.node_data_table = table( ...
        structGraph.node_labels(:), ...
        structGraph.full_node_indices, ...
        structGraph.node_radius, ...
        structGraph.node_thickness, ...
        'VariableNames', {'Label', 'FullNodeIndex', 'Radius', 'Thickness'});
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

function cyclePath = trace_pure_cycle(compNodes, G)
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
        cyclePath(end + 1) = next; %#ok<AGROW>
        prev = curr;
        curr = next;
        if curr == startNode
            break;
        end
    end
end

function edgeIds = path_to_edge_ids(path, edgeLookup)
    edgeIds = zeros(numel(path) - 1, 1);
    for i = 1:(numel(path) - 1)
        edgeIds(i) = edgeLookup(path(i), path(i + 1));
    end
    edgeIds = edgeIds(edgeIds > 0);
end

% =========================================================================
% Cleanup: contract short edges + merge close junctions
% =========================================================================

function [g, info] = cleanup_structural(g)
    info = struct('grid_spacing', NaN, 'min_edge_length', NaN, 'merge_radius', NaN, ...
        'iterations_used', 0, 'short_edge_contractions', 0, 'cluster_merges', 0);

    if g.num_nodes < 2 || isempty(g.edges_local)
        return;
    end

    h = estimate_characteristic_spacing(g.node_coords);
    if ~isfinite(h) || h <= 0
        return;
    end

    minEdgeLength = 1.5 * h;
    mergeRadius   = 1.5 * h;
    maxIterations = 4;
    info.grid_spacing = h;
    info.min_edge_length = minEdgeLength;
    info.merge_radius = mergeRadius;

    nShort = 0;
    nClust = 0;
    iterUsed = 0;
    for iter = 1:maxIterations
        iterUsed = iter;
        changed = false;

        while true
            pair = find_short_edge_pair(g, minEdgeLength);
            if isempty(pair)
                break;
            end
            g = merge_pair(g, pair);
            nShort = nShort + 1;
            changed = true;
        end

        while true
            pair = find_junction_cluster_pair(g, mergeRadius);
            if isempty(pair)
                break;
            end
            g = merge_pair(g, pair);
            nClust = nClust + 1;
            changed = true;
        end

        if ~changed
            break;
        end
    end

    info.iterations_used = iterUsed;
    info.short_edge_contractions = nShort;
    info.cluster_merges = nClust;
end

function h = estimate_characteristic_spacing(coords)
    xy = coords(:, 1:2);
    [xVals, ~] = cluster_axis_values(xy(:, 1));
    [yVals, ~] = cluster_axis_values(xy(:, 2));
    diffs = [diff(sort(xVals(:))); diff(sort(yVals(:)))];
    diffs = diffs(diffs > 0);
    if isempty(diffs)
        h = NaN;
    else
        h = median(diffs);
    end
end

function pair = find_short_edge_pair(g, minEdgeLength)
    pair = [];
    if isempty(g.edges_local)
        return;
    end
    [G, deg] = build_simple_graph(g);
    [~, order] = sort(g.edge_length(:), 'ascend');
    for k = 1:numel(order)
        idx = order(k);
        if g.edge_length(idx) >= minEdgeLength
            return;
        end
        u = g.edges_local(idx, 1);
        v = g.edges_local(idx, 2);
        if u == v || deg(u) == 1 || deg(v) == 1
            continue;
        end
        if deg(u) ~= 2 || deg(v) ~= 2 || share_common_neighbour(G, u, v)
            pair = [u, v];
            return;
        end
    end
end

function pair = find_junction_cluster_pair(g, mergeRadius)
    pair = [];
    if g.num_nodes < 2 || isempty(g.edges_local)
        return;
    end
    [G, deg] = build_simple_graph(g);
    compId = conncomp(G);
    cands = find(deg > 2);
    if numel(cands) < 2
        return;
    end
    xy = g.node_coords(:, 1:2);
    bestD = inf;
    for i = 1:(numel(cands) - 1)
        u = cands(i);
        for j = (i + 1):numel(cands)
            v = cands(j);
            if compId(u) ~= compId(v)
                continue;
            end
            d = norm(xy(u, :) - xy(v, :));
            if d <= mergeRadius && d < bestD
                bestD = d;
                pair = [u, v];
            end
        end
    end
end

function [G, deg] = build_simple_graph(g)
    if isempty(g.edges_local)
        G = graph([], [], [], g.num_nodes);
    else
        G = graph(g.edges_local(:, 1), g.edges_local(:, 2), [], g.num_nodes);
    end
    deg = degree(G);
end

function tf = share_common_neighbour(G, u, v)
    nbrU = neighbors(G, u);
    nbrV = neighbors(G, v);
    nbrU(nbrU == v) = [];
    nbrV(nbrV == u) = [];
    tf = ~isempty(intersect(nbrU, nbrV));
end

function g = merge_pair(g, pair)
    pair = unique(pair(:).', 'stable');
    pair = pair(pair >= 1 & pair <= g.num_nodes);
    if numel(pair) < 2
        return;
    end
    [~, deg] = build_simple_graph(g);
    rep = choose_representative(g, pair, deg);
    removeNodes = setdiff(pair, rep, 'stable');
    if isempty(removeNodes)
        return;
    end

    nE = size(g.edges_local, 1);
    for e = 1:nE
        src = g.edges_local(e, 1);
        dst = g.edges_local(e, 2);
        poly = g.edge_polyline{e};
        if ismember(src, removeNodes)
            g.edges_local(e, 1) = rep;
            if ~isempty(poly)
                poly(1, :) = g.node_coords(rep, :);
            end
        end
        if ismember(dst, removeNodes)
            g.edges_local(e, 2) = rep;
            if ~isempty(poly)
                poly(end, :) = g.node_coords(rep, :);
            end
        end
        g.edge_polyline{e} = poly;
    end

    removeMask = false(g.num_nodes, 1);
    removeMask(removeNodes) = true;
    g = compact_structural(g, removeMask);
end

function rep = choose_representative(g, pair, deg)
    degVals = deg(pair);
    cands = pair(degVals == max(degVals));
    if numel(cands) > 1
        thick = g.node_thickness(cands);
        thick(~isfinite(thick)) = -inf;
        cands = cands(thick == max(thick));
    end
    if numel(cands) == 1
        rep = cands(1);
        return;
    end
    xy = g.node_coords(cands, 1:2);
    totalDist = zeros(numel(cands), 1);
    for i = 1:numel(cands)
        delta = xy - xy(i, :);
        totalDist(i) = sum(sqrt(sum(delta .^ 2, 2)));
    end
    [~, bestRel] = min(totalDist);
    rep = cands(bestRel);
end

function g = compact_structural(g, removeMask)
    keepMask = ~removeMask(:);
    oldToNew = zeros(numel(keepMask), 1);
    oldToNew(keepMask) = 1:nnz(keepMask);

    g.node_coords    = g.node_coords(keepMask, :);
    g.node_labels    = g.node_labels(keepMask);
    g.num_nodes      = nnz(keepMask);
    g.full_node_indices = g.full_node_indices(keepMask);
    g.reduced_node_to_full_indices = g.full_node_indices;
    g.node_radius    = g.node_radius(keepMask);
    g.node_thickness = g.node_thickness(keepMask);
    g.node_data_table = g.node_data_table(keepMask, :);

    if isempty(g.edges_local)
        g.edge_index = zeros(2, 0);
        return;
    end
    g.edges_local      = oldToNew(g.edges_local);
    valid              = all(g.edges_local > 0, 2);
    g.edges_local      = g.edges_local(valid, :);
    g.edge_polyline    = g.edge_polyline(valid);
    g.edge_full_node_paths = g.edge_full_node_paths(valid);
    g.edge_length      = g.edge_length(valid);

    nonSelf = g.edges_local(:, 1) ~= g.edges_local(:, 2);
    g.edges_local      = g.edges_local(nonSelf, :);
    g.edge_polyline    = g.edge_polyline(nonSelf);
    g.edge_full_node_paths = g.edge_full_node_paths(nonSelf);
    g.edge_length      = g.edge_length(nonSelf);

    g = recompute_edge_lengths(g);
    g = deduplicate_edges(g);
    g.edge_index = g.edges_local.';
end

function g = recompute_edge_lengths(g)
    nE = size(g.edges_local, 1);
    for i = 1:nE
        poly = g.edge_polyline{i};
        if ~isempty(poly) && size(poly, 1) >= 2
            g.edge_length(i, 1) = path_length(poly);
        else
            d = g.node_coords(g.edges_local(i, 1), :) - g.node_coords(g.edges_local(i, 2), :);
            g.edge_length(i, 1) = norm(d);
        end
    end
end

function g = deduplicate_edges(g)
    if isempty(g.edges_local)
        return;
    end
    sortedPairs = sort(g.edges_local, 2);
    [~, ~, groupId] = unique(sortedPairs, 'rows', 'stable');
    keep = zeros(max(groupId), 1);
    for grp = 1:max(groupId)
        members = find(groupId == grp);
        if numel(members) == 1
            keep(grp) = members;
        else
            [~, bestRel] = max(g.edge_length(members));
            keep(grp) = members(bestRel(1));
        end
    end
    keep = sort(keep(keep > 0));
    g.edges_local      = g.edges_local(keep, :);
    g.edge_polyline    = g.edge_polyline(keep);
    g.edge_full_node_paths = g.edge_full_node_paths(keep);
    g.edge_length      = g.edge_length(keep);
end

% =========================================================================
% Edge thickness summary (min only)
% =========================================================================

function g = attach_edge_thickness_min(g, fullThickness)
    nE = size(g.edges_local, 1);
    edgeMin = nan(nE, 1);
    if nE == 0 || isempty(fullThickness)
        g.edge_thickness_min = edgeMin;
        return;
    end
    for i = 1:nE
        pathIdx = g.edge_full_node_paths{i};
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
    end
    g.edge_thickness_min = edgeMin;
end

% =========================================================================
% Misc
% =========================================================================

function len = path_length(coords)
    if size(coords, 1) < 2
        len = 0;
        return;
    end
    delta = diff(coords, 1, 1);
    len = sum(sqrt(sum(delta .^ 2, 2)));
end
