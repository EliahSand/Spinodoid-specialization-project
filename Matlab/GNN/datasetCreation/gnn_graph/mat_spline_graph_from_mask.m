function [matGraph, debugData] = mat_spline_graph_from_mask(occMask, xVals, yVals, varargin)
%MAT_SPLINE_GRAPH_FROM_MASK Build a compact MAT/spline control graph.
%
% [matGraph, debugData] = MAT_SPLINE_GRAPH_FROM_MASK(occMask, xVals, yVals, ...)
%
% Input:
%   occMask - logical occupancy image, rows correspond to yVals, columns to xVals
%   xVals   - physical x coordinate for each column
%   yVals   - physical y coordinate for each row
%
% Output graph schema:
%   node_features = [x, y, radius, degree, is_endpoint, is_junction, is_internal_control]
%
% The implementation follows the MAT idea used by Zhu et al.: medial points
% carry (x,y,radius), branches are traced as chains, and each branch is fitted
% by a compact cubic B-spline control polygon. Endpoint branches are pruned only
% when reconstructing the occupied mask from the remaining medial circles stays
% within the requested one-sided error tolerance.

    p = inputParser;
    p.addParameter('ErrorTolerance', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    p.addParameter('RidgeLow', 0.05, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('RidgeHigh', 0.45, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('LoGSigma', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('MinIslandPixels', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('MaxControlPointsPerBranch', 12, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    p.addParameter('MinControlPointsPerBranch', 4, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    p.addParameter('PruneEndpointBranches', true, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    occMask = logical(occMask);
    xVals = double(xVals(:));
    yVals = double(yVals(:));
    if size(occMask, 2) ~= numel(xVals) || size(occMask, 1) ~= numel(yVals)
        error('mat_spline_graph_from_mask:BadGrid', ...
            'occMask size must match numel(yVals) x numel(xVals).');
    end
    if ~any(occMask(:))
        error('mat_spline_graph_from_mask:EmptyMask', 'Occupancy mask is empty.');
    end

    h = estimate_grid_spacing(xVals, yVals);
    if ~isfinite(h) || h <= 0
        h = 1;
    end
    if isempty(opts.ErrorTolerance)
        opts.ErrorTolerance = 1.5 * h;
    end
    tolPix = opts.ErrorTolerance / h;

    distancePix = bwdist(~occMask, 'euclidean');
    radiusMap = distancePix * h;
    radiusMap(~occMask) = 0;

    ridgeness = compute_ridgeness(distancePix, occMask, opts.LoGSigma);
    maxBallSeedMask = local_distance_maxima(distancePix, occMask);
    strongRidgeMask = occMask & (ridgeness >= opts.RidgeHigh);
    ridgeGuideMask = occMask & (ridgeness >= opts.RidgeLow);

    componentAt = connected_components_binary(occMask, 8);
    skelMask = false(size(occMask));
    nComp = max(componentAt(:));
    for compId = 1:nComp
        compMask = componentAt == compId;
        compSkel = skeletonize_component(compMask);
        if nnz(compSkel) < opts.MinIslandPixels
            compIdx = find(compMask);
            [~, pickRel] = max(distancePix(compIdx));
            compSkel(compIdx(pickRel)) = true;
        end

        % Keep ridge/maximal-ball pixels that already lie on the thin medial
        % skeleton. They do not thicken the skeleton, but they survive pruning
        % as reconstruction-important medial samples.
        compSeeds = compSkel & (maxBallSeedMask | strongRidgeMask | ridgeGuideMask);
        if any(compSeeds(:))
            compSkel = compSkel | compSeeds;
            compSkel = skeletonize_component(compSkel);
        end
        skelMask = skelMask | compSkel;
    end

    if logical(opts.PruneEndpointBranches)
        skelMask = prune_endpoint_branches(skelMask, occMask, distancePix, tolPix);
    end

    reconMask = reconstruct_from_medial_circles(skelMask, distancePix);
    finalErrorPix = one_sided_missing_error(occMask, reconMask);

    skelGraph = skeleton_graph_from_mask(skelMask, componentAt, radiusMap, xVals, yVals);
    branches = trace_skeleton_branches(skelGraph);
    matGraph = build_control_graph_from_branches(skelGraph, branches, opts, h);

    matGraph.kind = 'mat_spline_control_graph';
    matGraph.schema_version = 4;
    matGraph.feature_names = {'x', 'y', 'radius', 'degree', ...
        'is_endpoint', 'is_junction', 'is_internal_control'};
    matGraph.error_tolerance = opts.ErrorTolerance;
    matGraph.grid_spacing = h;
    matGraph.reconstruction_error = finalErrorPix * h;
    matGraph.skeleton_num_nodes = skelGraph.num_nodes;
    matGraph.skeleton_num_edges = size(skelGraph.edges_local, 1);
    matGraph.ridge_low = opts.RidgeLow;
    matGraph.ridge_high = opts.RidgeHigh;

    debugData = struct();
    debugData.occupancy_mask = occMask;
    debugData.radius_map = radiusMap;
    debugData.distance_pixels = distancePix;
    debugData.ridgeness = ridgeness;
    debugData.max_ball_seed_mask = maxBallSeedMask;
    debugData.strong_ridge_mask = strongRidgeMask;
    debugData.ridge_guide_mask = ridgeGuideMask;
    debugData.component_at_grid = componentAt;
    debugData.skeleton_mask = skelMask;
    debugData.reconstruction_mask = reconMask;
    debugData.reconstruction_error = matGraph.reconstruction_error;
    debugData.error_tolerance = opts.ErrorTolerance;
    debugData.occupancy_component_count = nComp;
    debugData.skeleton_component_count = max(connected_components_binary(skelMask, 8), [], 'all');
    debugData.skeleton_graph = skelGraph;
end

function h = estimate_grid_spacing(xVals, yVals)
    dx = diff(sort(xVals(:)));
    dy = diff(sort(yVals(:)));
    vals = [dx(dx > 0); dy(dy > 0)];
    if isempty(vals)
        h = NaN;
    else
        h = median(vals);
    end
end

function rdg = compute_ridgeness(distancePix, mask, sigma)
    filterSize = max(3, 2 * ceil(3 * sigma) + 1);
    logKernel = fspecial('log', filterSize, sigma);
    rdg = -imfilter(distancePix, logKernel, 'replicate', 'same');
    rdg(~mask) = 0;
    vals = rdg(mask);
    vals = vals(isfinite(vals));
    if isempty(vals)
        rdg(:) = 0;
        return;
    end
    lo = min(vals);
    hi = max(vals);
    rdg = (rdg - lo) / max(hi - lo, eps);
    rdg(~mask) = 0;
end

function maxMask = local_distance_maxima(distancePix, mask)
    padded = -inf(size(distancePix) + 2);
    padded(2:end-1, 2:end-1) = distancePix;
    maxNbr = -inf(size(distancePix));
    for dr = -1:1
        for dc = -1:1
            if dr == 0 && dc == 0
                continue;
            end
            block = padded((2 + dr):(end - 1 + dr), (2 + dc):(end - 1 + dc));
            maxNbr = max(maxNbr, block);
        end
    end
    maxMask = mask & distancePix >= maxNbr & distancePix > 0;
end

function skel = skeletonize_component(mask)
    mask = logical(mask);
    if ~any(mask(:))
        skel = false(size(mask));
        return;
    end
    try
        skel = bwskel(mask);
    catch
        skel = bwmorph(mask, 'skel', Inf);
    end
    skel = logical(skel);
end

function labels = connected_components_binary(mask, conn)
    if nargin < 2
        conn = 8;
    end
    mask = logical(mask);
    labels = zeros(size(mask));
    if ~any(mask(:))
        return;
    end
    if conn == 4
        neigh = [-1 0; 1 0; 0 -1; 0 1];
    else
        neigh = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];
    end
    [rows, cols] = find(mask);
    queue = zeros(nnz(mask), 2);
    current = 0;
    for seed = 1:numel(rows)
        r0 = rows(seed);
        c0 = cols(seed);
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
            for k = 1:size(neigh, 1)
                rr = r + neigh(k, 1);
                cc = c + neigh(k, 2);
                if rr < 1 || rr > size(mask, 1) || cc < 1 || cc > size(mask, 2)
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

function pruned = prune_endpoint_branches(skelMask, occMask, distancePix, tolPix)
    pruned = logical(skelMask);
    maxIterations = max(1, nnz(pruned));
    for iter = 1:maxIterations
        skelGraph = skeleton_graph_from_mask(pruned, zeros(size(pruned)), distancePix, ...
            1:size(pruned, 2), 1:size(pruned, 1));
        if skelGraph.num_nodes == 0
            break;
        end
        deg = degree(skelGraph.graph);
        endpoints = find(deg == 1);
        changed = false;
        for e = endpoints(:).'
            path = trace_endpoint_path(e, skelGraph.graph, deg);
            if numel(path) < 2
                continue;
            end
            if deg(path(end)) == 1
                removePath = path(1:end-1);
            else
                removePath = path(1:end-1);
            end
            if isempty(removePath)
                continue;
            end
            candidate = pruned;
            candidate(skelGraph.grid_linear_index(removePath)) = false;
            if ~any(candidate(:))
                continue;
            end
            recon = reconstruct_from_medial_circles(candidate, distancePix);
            errPix = one_sided_missing_error(occMask, recon);
            if errPix <= tolPix + eps
                pruned = candidate;
                changed = true;
                break;
            end
        end
        if ~changed
            break;
        end
    end
end

function path = trace_endpoint_path(endpoint, G, deg)
    path = endpoint;
    nbr = neighbors(G, endpoint);
    if isempty(nbr)
        return;
    end
    prev = endpoint;
    curr = nbr(1);
    path(end + 1) = curr;
    while deg(curr) == 2
        nbrs = neighbors(G, curr).';
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            break;
        end
        next = nbrs(1);
        path(end + 1) = next; %#ok<AGROW>
        prev = curr;
        curr = next;
    end
end

function recon = reconstruct_from_medial_circles(skelMask, distancePix)
    recon = false(size(skelMask));
    [rows, cols] = find(skelMask);
    [ccGrid, rrGrid] = meshgrid(1:size(skelMask, 2), 1:size(skelMask, 1));
    for k = 1:numel(rows)
        r = rows(k);
        c = cols(k);
        rad = distancePix(r, c);
        if rad <= 0
            recon(r, c) = true;
            continue;
        end
        rLo = max(1, floor(r - rad));
        rHi = min(size(skelMask, 1), ceil(r + rad));
        cLo = max(1, floor(c - rad));
        cHi = min(size(skelMask, 2), ceil(c + rad));
        local = (rrGrid(rLo:rHi, cLo:cHi) - r) .^ 2 + ...
            (ccGrid(rLo:rHi, cLo:cHi) - c) .^ 2 <= rad ^ 2 + eps;
        recon(rLo:rHi, cLo:cHi) = recon(rLo:rHi, cLo:cHi) | local;
    end
end

function errPix = one_sided_missing_error(occMask, reconMask)
    missing = occMask & ~reconMask;
    if ~any(missing(:))
        errPix = 0;
        return;
    end
    if ~any(reconMask(:))
        errPix = inf;
        return;
    end
    distToRecon = bwdist(reconMask, 'euclidean');
    errPix = max(distToRecon(missing));
end

function skelGraph = skeleton_graph_from_mask(skelMask, componentAt, radiusMap, xVals, yVals)
    skelMask = logical(skelMask);
    skelLinear = find(skelMask);
    nSkel = numel(skelLinear);
    skelLocalAt = zeros(size(skelMask));
    skelLocalAt(skelLinear) = 1:nSkel;

    coords = zeros(nSkel, 2);
    radii = zeros(nSkel, 1);
    for i = 1:nSkel
        [r, c] = ind2sub(size(skelMask), skelLinear(i));
        coords(i, :) = [xVals(c), yVals(r)];
        radii(i) = radiusMap(r, c);
    end

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
            j = skelLocalAt(rr, cc);
            if j == 0
                continue;
            end
            if any(componentAt(:)) && componentAt(r, c) ~= 0 && componentAt(rr, cc) ~= componentAt(r, c)
                continue;
            end
            if abs(dirs(d, 1)) == 1 && abs(dirs(d, 2)) == 1 && any(componentAt(:))
                compId = componentAt(r, c);
                if compId ~= 0 && componentAt(r, cc) ~= compId && componentAt(rr, c) ~= compId
                    continue;
                end
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

    skelGraph = struct();
    skelGraph.kind = 'mat_skeleton_graph';
    skelGraph.num_nodes = nSkel;
    skelGraph.node_coords = coords;
    skelGraph.node_radius = radii;
    skelGraph.edges_local = edges;
    skelGraph.edge_index = edges.';
    skelGraph.edge_length = lengths;
    skelGraph.graph = G;
    skelGraph.grid_linear_index = skelLinear(:);
end

function branches = trace_skeleton_branches(skelGraph)
    branches = struct('path', {}, 'is_cycle', {});
    n = skelGraph.num_nodes;
    if n == 0
        return;
    end
    G = skelGraph.graph;
    deg = degree(G);
    edges = skelGraph.edges_local;
    if isempty(edges)
        for i = 1:n
            branches(end + 1).path = i; %#ok<AGROW>
            branches(end).is_cycle = false;
        end
        return;
    end
    nEdge = size(edges, 1);
    edgeLookup = sparse([edges(:, 1); edges(:, 2)], [edges(:, 2); edges(:, 1)], ...
        [(1:nEdge).'; (1:nEdge).'], n, n);
    visited = false(nEdge, 1);
    compId = conncomp(G);

    for c = 1:max(compId)
        compNodes = find(compId == c);
        keyNodes = compNodes(deg(compNodes) ~= 2);
        if isempty(keyNodes)
            path = trace_pure_cycle(compNodes, G, edgeLookup);
            if numel(path) >= 2
                branches(end + 1).path = path; %#ok<AGROW>
                branches(end).is_cycle = true;
                visited(path_to_edge_ids(path, edgeLookup)) = true;
            end
            continue;
        end
        for u = sort(keyNodes(:)).'
            nbrs = sort(neighbors(G, u)).';
            if isempty(nbrs)
                branches(end + 1).path = u; %#ok<AGROW>
                branches(end).is_cycle = false;
                continue;
            end
            for v = nbrs
                eIdx = edgeLookup(u, v);
                if eIdx == 0 || visited(eIdx)
                    continue;
                end
                path = trace_chain(u, v, G, deg, edgeLookup, visited);
                edgeIds = path_to_edge_ids(path, edgeLookup);
                visited(edgeIds) = true;
                branches(end + 1).path = path; %#ok<AGROW>
                branches(end).is_cycle = path(1) == path(end);
            end
        end
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
        nxt = nbrs(1);
        eIdx = edgeLookup(curr, nxt);
        if eIdx == 0 || (visited(eIdx) && nxt ~= startNode)
            break;
        end
        path(end + 1) = nxt; %#ok<AGROW>
        prev = curr;
        curr = nxt;
    end
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
    maxSteps = numel(compNodes) + 2;
    while numel(cyclePath) <= maxSteps
        nbrs = sort(neighbors(G, curr)).';
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            break;
        end
        nxt = nbrs(1);
        if edgeLookup(curr, nxt) == 0
            break;
        end
        cyclePath(end + 1) = nxt; %#ok<AGROW>
        prev = curr;
        curr = nxt;
        if curr == startNode
            break;
        end
    end
end

function edgeIds = path_to_edge_ids(path, edgeLookup)
    edgeIds = zeros(max(0, numel(path) - 1), 1);
    for i = 1:(numel(path) - 1)
        edgeIds(i) = edgeLookup(path(i), path(i + 1));
    end
    edgeIds = edgeIds(edgeIds > 0);
end

function matGraph = build_control_graph_from_branches(skelGraph, branches, opts, h)
    coords = zeros(0, 2);
    radius = zeros(0, 1);
    keyRole = zeros(0, 1); % 0 internal control, 1 endpoint, 2 junction
    reducedNodeForSkel = zeros(skelGraph.num_nodes, 1);
    edges = zeros(0, 2);
    edgeSet = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    Gskel = skelGraph.graph;
    skelDegree = degree(Gskel);
    branchMeta = struct('control_indices', {}, 'source_skeleton_indices', {}, ...
        'knots', {}, 'degree', {}, 'fit_error', {}, 'path_length', {}, ...
        'radius_min', {}, 'radius_mean', {}, 'radius_max', {}, 'radius_std', {}, ...
        'is_cycle', {});

    for b = 1:numel(branches)
        path = branches(b).path(:).';
        if isempty(path)
            continue;
        end
        if numel(path) == 1
            nodeId = ensure_key_node(path(1));
            add_branch_meta(nodeId, path, [], 0, 0, 0, false);
            continue;
        end

        points3 = [skelGraph.node_coords(path, :), skelGraph.node_radius(path)];
        fit = fit_branch_spline(points3, opts.ErrorTolerance, ...
            opts.MinControlPointsPerBranch, opts.MaxControlPointsPerBranch);

        if branches(b).is_cycle
            ctrl = fit.control_points;
            if size(ctrl, 1) > 3 && norm(ctrl(end, 1:2) - ctrl(1, 1:2)) <= max(h, eps)
                ctrl = ctrl(1:end-1, :);
            end
            ctrlIds = zeros(1, size(ctrl, 1));
            for i = 1:size(ctrl, 1)
                ctrlIds(i) = add_control_node(ctrl(i, 1:2), ctrl(i, 3), 0);
            end
            add_edges_for_branch(ctrlIds, true, points3);
        else
            src = ensure_key_node(path(1));
            dst = ensure_key_node(path(end));
            ctrlIds = src;
            ctrl = fit.control_points;
            for i = 2:(size(ctrl, 1) - 1)
                ctrlIds(end + 1) = add_control_node(ctrl(i, 1:2), ctrl(i, 3), 0); %#ok<AGROW>
            end
            ctrlIds(end + 1) = dst;
            ctrlIds = ensure_nonduplicate_direct_edge(ctrlIds, points3);
            add_edges_for_branch(ctrlIds, false, points3);
        end

        add_branch_meta(ctrlIds, path, fit.knots, fit.degree, fit.fit_error, ...
            path_length(points3(:, 1:2)), branches(b).is_cycle);
    end

    n = size(coords, 1);
    if isempty(edges)
        G = graph([], [], [], n);
    else
        G = graph(edges(:, 1), edges(:, 2), [], n);
    end
    deg = degree(G);
    isEndpoint = keyRole == 1 | deg == 1;
    isJunction = keyRole == 2 | deg > 2;
    isInternal = keyRole == 0;
    nodeFeatures = [coords(:, 1:2), radius(:), double(deg(:)), ...
        double(isEndpoint(:)), double(isJunction(:)), double(isInternal(:))];

    matGraph = struct();
    matGraph.num_nodes = n;
    matGraph.num_edges = size(edges, 1);
    matGraph.node_coords = coords;
    matGraph.node_radius = radius(:);
    matGraph.node_degree = deg(:);
    matGraph.node_features = nodeFeatures;
    matGraph.edges_local = edges;
    matGraph.edge_index = int32(edges.');
    matGraph.graph = G;
    matGraph.branches = branchMeta;
    matGraph.branch_count = numel(branchMeta);

    function nodeId = ensure_key_node(skelIdx)
        nodeId = reducedNodeForSkel(skelIdx);
        if nodeId > 0
            return;
        end
        if skelDegree(skelIdx) == 1
            role = 1;
        elseif skelDegree(skelIdx) > 2
            role = 2;
        else
            role = 0;
        end
        nodeId = add_control_node(skelGraph.node_coords(skelIdx, :), skelGraph.node_radius(skelIdx), role);
        reducedNodeForSkel(skelIdx) = nodeId;
    end

    function nodeId = add_control_node(xy, r, role)
        nodeId = size(coords, 1) + 1;
        coords(nodeId, :) = xy; %#ok<AGROW>
        radius(nodeId, 1) = max(0, r); %#ok<AGROW>
        keyRole(nodeId, 1) = role; %#ok<AGROW>
    end

    function ctrlIds = ensure_nonduplicate_direct_edge(ctrlIds, points3)
        if numel(ctrlIds) ~= 2
            return;
        end
        if ~has_edge(ctrlIds(1), ctrlIds(2))
            return;
        end
        midIdx = max(1, round(size(points3, 1) / 2));
        midNode = add_control_node(points3(midIdx, 1:2), points3(midIdx, 3), 0);
        ctrlIds = [ctrlIds(1), midNode, ctrlIds(2)];
    end

    function add_edges_for_branch(ctrlIds, closeLoop, points3)
        if numel(ctrlIds) < 2
            return;
        end
        for i = 1:(numel(ctrlIds) - 1)
            add_unique_edge(ctrlIds(i), ctrlIds(i + 1), points3);
        end
        if closeLoop && numel(ctrlIds) > 2
            add_unique_edge(ctrlIds(end), ctrlIds(1), points3);
        end
    end

    function add_unique_edge(u, v, points3)
        if u == v
            return;
        end
        key = edge_key(u, v);
        if isKey(edgeSet, key)
            midIdx = max(1, round(size(points3, 1) / 2));
            midNode = add_control_node(points3(midIdx, 1:2), points3(midIdx, 3), 0);
            add_unique_edge(u, midNode, points3);
            add_unique_edge(midNode, v, points3);
            return;
        end
        edgeSet(key) = true;
        edges(end + 1, :) = [u, v]; %#ok<AGROW>
    end

    function tf = has_edge(u, v)
        tf = isKey(edgeSet, edge_key(u, v));
    end

    function add_branch_meta(ctrlIds, path, knots, degree, fitError, len, isCycle)
        rr = skelGraph.node_radius(path);
        branchMeta(end + 1).control_indices = int32(ctrlIds(:)); %#ok<AGROW>
        branchMeta(end).source_skeleton_indices = int32(path(:));
        branchMeta(end).knots = knots(:).';
        branchMeta(end).degree = degree;
        branchMeta(end).fit_error = fitError;
        branchMeta(end).path_length = len;
        branchMeta(end).radius_min = min(rr);
        branchMeta(end).radius_mean = mean(rr);
        branchMeta(end).radius_max = max(rr);
        branchMeta(end).radius_std = std(rr);
        branchMeta(end).is_cycle = logical(isCycle);
    end
end

function key = edge_key(u, v)
    if u > v
        tmp = u;
        u = v;
        v = tmp;
    end
    key = sprintf('%d_%d', u, v);
end

function fit = fit_branch_spline(points3, errorTol, minCtrl, maxCtrl)
    nData = size(points3, 1);
    if nData <= 1
        fit = struct('control_points', points3, 'knots', [0 1], ...
            'degree', 0, 'fit_error', 0);
        return;
    end

    t = chord_parameters(points3(:, 1:2));
    maxCtrl = min(maxCtrl, nData);
    minCtrl = min(max(2, minCtrl), maxCtrl);
    best = [];
    for nCtrl = minCtrl:maxCtrl
        degree = min(3, nCtrl - 1);
        knots = open_uniform_knots(nCtrl, degree);
        B = bspline_basis_matrix(t, nCtrl, degree, knots);
        C = fit_control_points(B, points3);
        approx = B * C;
        err = max(sqrt(sum((points3 - approx) .^ 2, 2)));
        best = struct('control_points', C, 'knots', knots, ...
            'degree', degree, 'fit_error', err);
        if err <= errorTol + eps
            break;
        end
    end
    fit = best;
end

function t = chord_parameters(xy)
    if size(xy, 1) == 1
        t = 0;
        return;
    end
    ds = sqrt(sum(diff(xy, 1, 1) .^ 2, 2));
    s = [0; cumsum(ds)];
    if s(end) <= eps
        t = linspace(0, 1, size(xy, 1)).';
    else
        t = s / s(end);
    end
    t(1) = 0;
    t(end) = 1;
end

function knots = open_uniform_knots(nCtrl, degree)
    nInterior = nCtrl - degree - 1;
    if nInterior > 0
        interior = (1:nInterior) / (nInterior + 1);
    else
        interior = zeros(1, 0);
    end
    knots = [zeros(1, degree + 1), interior, ones(1, degree + 1)];
end

function B = bspline_basis_matrix(t, nCtrl, degree, knots)
    t = t(:);
    B = zeros(numel(t), nCtrl);
    for i = 1:nCtrl
        B(:, i) = bspline_basis(i, degree, t, knots, nCtrl);
    end
    B(t == 1, :) = 0;
    B(t == 1, nCtrl) = 1;
end

function N = bspline_basis(i, degree, t, knots, nCtrl)
    if degree == 0
        N = double(knots(i) <= t & t < knots(i + 1));
        if i == nCtrl
            N(t == 1) = 1;
        end
        return;
    end
    denomA = knots(i + degree) - knots(i);
    denomB = knots(i + degree + 1) - knots(i + 1);
    A = zeros(size(t));
    B = zeros(size(t));
    if denomA > eps
        A = ((t - knots(i)) / denomA) .* bspline_basis(i, degree - 1, t, knots, nCtrl);
    end
    if denomB > eps
        B = ((knots(i + degree + 1) - t) / denomB) .* bspline_basis(i + 1, degree - 1, t, knots, nCtrl);
    end
    N = A + B;
end

function C = fit_control_points(B, points3)
    nCtrl = size(B, 2);
    if nCtrl <= 2
        C = B \ points3;
        return;
    end
    C = zeros(nCtrl, size(points3, 2));
    C(1, :) = points3(1, :);
    C(end, :) = points3(end, :);
    rhs = points3 - B(:, [1, nCtrl]) * C([1, nCtrl], :);
    C(2:end-1, :) = B(:, 2:end-1) \ rhs;
end

function len = path_length(xy)
    if size(xy, 1) < 2
        len = 0;
        return;
    end
    len = sum(sqrt(sum(diff(xy, 1, 1) .^ 2, 2)));
end
