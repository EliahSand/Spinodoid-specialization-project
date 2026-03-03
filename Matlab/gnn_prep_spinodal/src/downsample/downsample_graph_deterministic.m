function downGraph = downsample_graph_deterministic(fullGraph, varargin)
%DOWNSAMPLE_GRAPH_DETERMINISTIC Deterministic graph downsampling for GNN prep.
%
% downGraph = DOWNSAMPLE_GRAPH_DETERMINISTIC(fullGraph, ...)
%
% Name-Value options:
%   'DetailLevel'         : [0,1], fraction of interior nodes to keep (default 0.35)
%   'TargetNumNodes'      : explicit target node count (optional)
%   'PreserveConnectivity': if true, add shortest-path connector nodes (default true)
%   'NodeSelectionMethod' : 'farthest'|'grid'|'hybrid' (default 'farthest')
%   'HybridAnchorFraction': [0,1], anchor ratio in hybrid mode (default 0.35)
%   'OutputMode'          : 'mesh' or 'line' (default 'mesh')
%   'LineNodeBudget'      : max nodes kept in line mode (optional)
%   'BoundaryKeepRatio'   : [0,1], boundary-node keep ratio in mesh mode (default 1)
%   'BoundaryStraightTolDeg' : corner detection tolerance in degrees (default 8)
%
% Guarantees:
% - Deterministic/stable for same inputs (no randomness).
% - Boundary retention in mesh mode is controlled by 'BoundaryKeepRatio':
%   1.0 keeps all boundary nodes; values < 1 simplify boundary chains.
% - Simplified boundary chains are re-linked with shortcut edges so the
%   boundary stays connected after node removal.

    p = inputParser;
    p.addParameter('DetailLevel', 0.35, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('TargetNumNodes', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    p.addParameter('PreserveConnectivity', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('NodeSelectionMethod', 'farthest', @(x) ischar(x) || isstring(x));
    p.addParameter('HybridAnchorFraction', 0.35, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('OutputMode', 'mesh', @(x) ischar(x) || isstring(x));
    p.addParameter('LineNodeBudget', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    p.addParameter('BoundaryKeepRatio', 1, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.addParameter('BoundaryStraightTolDeg', 8, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 180);
    p.parse(varargin{:});
    opts = p.Results;
    nodeMethod = lower(char(string(opts.NodeSelectionMethod)));
    nodeMethod = validatestring(nodeMethod, {'farthest', 'grid', 'hybrid'});
    outputMode = lower(char(string(opts.OutputMode)));
    outputMode = validatestring(outputMode, {'mesh', 'line'});

    required = {'node_labels', 'node_coords', 'edges_local', 'boundary_mask'};
    for i = 1:numel(required)
        if ~isfield(fullGraph, required{i})
            error('downsample_graph_deterministic:MissingField', ...
                'fullGraph missing required field: %s', required{i});
        end
    end

    nodeLabels = fullGraph.node_labels(:);
    nodeCoords = fullGraph.node_coords;
    nNodes = numel(nodeLabels);

    boundaryMask = logical(fullGraph.boundary_mask(:));
    if numel(boundaryMask) ~= nNodes
        error('downsample_graph_deterministic:BadBoundaryMask', ...
            'boundary_mask length does not match node count.');
    end

    edgesLocal = fullGraph.edges_local;
    if isfield(fullGraph, 'graph') && isa(fullGraph.graph, 'graph')
        G = fullGraph.graph;
    else
        if isempty(edgesLocal)
            G = graph([], [], [], nNodes);
        else
            G = graph(edgesLocal(:, 1), edgesLocal(:, 2), [], nNodes);
        end
    end

    boundaryIdx = find(boundaryMask);
    interiorIdx = find(~boundaryMask);
    nBoundary = numel(boundaryIdx);
    nInterior = numel(interiorIdx);

    selectedBoundary = boundaryIdx(:);
    selectedBoundaryEdgesFull = zeros(0, 2);
    if opts.BoundaryKeepRatio < 1 && ~isempty(selectedBoundary)
        boundaryEdgesFull = get_boundary_edges(fullGraph, boundaryMask);
        if ~isempty(boundaryEdgesFull)
            [simplifiedBoundary, simplifiedBoundaryEdges] = simplify_polyline_components( ...
                boundaryEdgesFull, nodeCoords, nodeLabels, ...
                opts.BoundaryKeepRatio, opts.BoundaryStraightTolDeg);
            simplifiedBoundary = intersect(simplifiedBoundary(:), boundaryIdx(:), 'stable');
            if ~isempty(simplifiedBoundary)
                selectedBoundary = simplifiedBoundary;
                if ~isempty(simplifiedBoundaryEdges)
                    selectedBoundaryEdgesFull = simplifiedBoundaryEdges;
                end
            end
        end
    end

    nBoundarySelected = numel(selectedBoundary);
    if isempty(opts.TargetNumNodes)
        targetNum = nBoundarySelected + round(opts.DetailLevel * nInterior);
    else
        targetNum = round(opts.TargetNumNodes);
    end
    targetNum = max(targetNum, nBoundarySelected);
    targetNum = min(targetNum, nNodes);

    selected = selectedBoundary(:);
    if isempty(selected) && nNodes > 0
        [~, seedPos] = min(nodeLabels);
        selected = seedPos;
    end

    nInteriorTarget = max(targetNum - numel(selected), 0);

    pickedInterior = select_interior_nodes( ...
        G, interiorIdx(:), selected(:), nInteriorTarget, ...
        nodeCoords, nodeLabels, nodeMethod, opts.HybridAnchorFraction);

    selectedMesh = unique([selected(:); pickedInterior(:)], 'sorted');
    selectedBeforeConnectivity = selectedMesh;

    if logical(opts.PreserveConnectivity)
        selectedMesh = augment_connectivity_with_shortest_paths(G, selectedMesh, nodeLabels);
    end
    selectedMesh = unique(selectedMesh(:), 'sorted');

    mapOldToNew = zeros(nNodes, 1);
    mapOldToNew(selectedMesh) = 1:numel(selectedMesh);

    edgesMesh = build_subset_edges(edgesLocal, selectedMesh, mapOldToNew);
    nBoundaryShortcutEdges = 0;
    if ~isempty(selectedBoundaryEdgesFull)
        % Preserve a continuous boundary line after dropping intermediate boundary nodes.
        boundaryShortcutEdges = build_subset_edges(selectedBoundaryEdgesFull, selectedMesh, mapOldToNew);
        if ~isempty(boundaryShortcutEdges)
            nBoundaryShortcutEdges = size(boundaryShortcutEdges, 1);
            edgesMesh = [edgesMesh; boundaryShortcutEdges]; %#ok<AGROW>
            edgesMesh = sort(edgesMesh, 2);
            edgesMesh = unique(edgesMesh, 'rows', 'stable');
        end
    end

    if isempty(edgesMesh)
        Gmesh = graph([], [], [], numel(selectedMesh));
    else
        Gmesh = graph(edgesMesh(:, 1), edgesMesh(:, 2), [], numel(selectedMesh));
    end

    selected = selectedMesh;
    edgesNew = edgesMesh;
    Gdown = Gmesh;
    if strcmp(outputMode, 'line')
        boundaryEdgesMesh = map_boundary_edges_to_subset(fullGraph, selectedMesh, nNodes);
        if isempty(boundaryEdgesMesh)
            % Fallback: simplify all current edges component-wise.
            [lineNodesMesh, lineEdgesMesh] = simplify_polyline_components( ...
                edgesMesh, nodeCoords(selectedMesh, :), nodeLabels(selectedMesh), ...
                derive_line_keep_ratio(opts.DetailLevel, opts.LineNodeBudget, numel(selectedMesh)), ...
                opts.BoundaryStraightTolDeg);
        else
            [lineNodesMesh, lineEdgesMesh] = simplify_polyline_components( ...
                boundaryEdgesMesh, nodeCoords(selectedMesh, :), nodeLabels(selectedMesh), ...
                derive_line_keep_ratio(opts.DetailLevel, opts.LineNodeBudget, ...
                    numel(unique(boundaryEdgesMesh(:)))), ...
                opts.BoundaryStraightTolDeg);
        end

        if isempty(lineNodesMesh)
            lineNodesMesh = 1:numel(selectedMesh);
            lineEdgesMesh = edgesMesh;
        end

        selected = selectedMesh(lineNodesMesh(:));
        mapLine = zeros(numel(selectedMesh), 1);
        mapLine(lineNodesMesh(:)) = 1:numel(lineNodesMesh);
        if isempty(lineEdgesMesh)
            edgesNew = zeros(0, 2);
        else
            edgesNew = [mapLine(lineEdgesMesh(:, 1)), mapLine(lineEdgesMesh(:, 2))];
            edgesNew = sort(edgesNew, 2);
            edgesNew = unique(edgesNew(edgesNew(:, 1) > 0 & edgesNew(:, 2) > 0, :), 'rows', 'stable');
        end
        if isempty(edgesNew)
            Gdown = graph([], [], [], numel(selected));
        else
            Gdown = graph(edgesNew(:, 1), edgesNew(:, 2), [], numel(selected));
        end
    end

    downGraph = struct();
    downGraph.kind = 'downsampled_graph';
    downGraph.parent_kind = get_default(fullGraph, 'kind', 'full_reference_graph');
    downGraph.detail_level = opts.DetailLevel;
    downGraph.target_num_nodes = targetNum;
    downGraph.actual_num_nodes = numel(selected);
    downGraph.preserve_connectivity = logical(opts.PreserveConnectivity);
    downGraph.node_selection_method = nodeMethod;
    downGraph.hybrid_anchor_fraction = opts.HybridAnchorFraction;
    downGraph.boundary_keep_ratio = opts.BoundaryKeepRatio;
    downGraph.boundary_straight_tol_deg = opts.BoundaryStraightTolDeg;
    downGraph.output_mode = outputMode;
    downGraph.line_node_budget = opts.LineNodeBudget;
    downGraph.is_line_representation = strcmp(outputMode, 'line');

    downGraph.selected_full_indices = selected;
    downGraph.selected_full_indices_pre_connectivity = selectedBeforeConnectivity;
    downGraph.selected_full_indices_mesh = selectedMesh;
    downGraph.selected_full_labels = nodeLabels(selected);

    downGraph.node_labels = nodeLabels(selected);
    downGraph.node_coords = nodeCoords(selected, :);
    downGraph.num_nodes = numel(selected);

    downGraph.edges_local = edgesNew;
    downGraph.graph = Gdown;

    downGraph.boundary_mask = boundaryMask(selected);
    downGraph.boundary_labels = downGraph.node_labels(downGraph.boundary_mask);
    if isfield(fullGraph, 'boundary_labels')
        downGraph.preserved_all_boundary_nodes = all(ismember(fullGraph.boundary_labels, downGraph.node_labels));
    else
        downGraph.preserved_all_boundary_nodes = all(ismember(nodeLabels(boundaryIdx), downGraph.node_labels));
    end

    downGraph.downsample_info = struct( ...
        'n_full_nodes', nNodes, ...
        'n_full_boundary_nodes', nBoundary, ...
        'n_boundary_nodes_after_boundary_downsample', nBoundarySelected, ...
        'n_full_interior_nodes', nInterior, ...
        'node_selection_method', nodeMethod, ...
        'hybrid_anchor_fraction', opts.HybridAnchorFraction, ...
        'n_boundary_shortcut_edges', nBoundaryShortcutEdges, ...
        'n_initial_selected_nodes', numel(selectedBeforeConnectivity), ...
        'n_added_for_connectivity', numel(selectedMesh) - numel(selectedBeforeConnectivity), ...
        'n_mesh_selected_nodes', numel(selectedMesh), ...
        'n_final_selected_nodes', numel(selected), ...
        'n_removed_by_line_projection', numel(selectedMesh) - numel(selected));

    if isfield(fullGraph, 'node_data_table') && istable(fullGraph.node_data_table) ...
            && height(fullGraph.node_data_table) == nNodes
        downGraph.node_data_table = fullGraph.node_data_table(selected, :);
        if isfield(fullGraph, 'node_attr_names')
            downGraph.node_attr_names = fullGraph.node_attr_names;
        end
    end
end

function picked = select_interior_nodes(G, interiorIdx, seedIdx, nToPick, nodeCoords, nodeLabels, method, hybridAnchorFraction)
    picked = zeros(0, 1);
    if nToPick <= 0 || isempty(interiorIdx)
        return;
    end

    method = validatestring(lower(char(string(method))), {'farthest', 'grid', 'hybrid'});
    switch method
        case 'grid'
            picked = grid_interior_sampling(G, interiorIdx, seedIdx, nToPick, nodeCoords, nodeLabels);
        case 'hybrid'
            picked = hybrid_interior_sampling( ...
                G, interiorIdx, seedIdx, nToPick, nodeCoords, nodeLabels, hybridAnchorFraction);
        otherwise
            picked = farthest_interior_sampling(G, interiorIdx, seedIdx, nToPick, nodeLabels);
    end

    picked = unique(picked(:), 'stable');
    if numel(picked) > nToPick
        picked = picked(1:nToPick);
        return;
    end

    if numel(picked) < nToPick
        missing = nToPick - numel(picked);
        remaining = setdiff(unique(interiorIdx(:), 'sorted'), picked, 'stable');
        if ~isempty(remaining)
            extra = farthest_interior_sampling(G, remaining, [seedIdx(:); picked(:)], missing, nodeLabels);
            picked = unique([picked(:); extra(:)], 'stable');
        end
    end
end

function picked = hybrid_interior_sampling(G, interiorIdx, seedIdx, nToPick, nodeCoords, nodeLabels, anchorFraction)
    picked = zeros(0, 1);
    if nToPick <= 0 || isempty(interiorIdx)
        return;
    end

    candidates = unique(interiorIdx(:), 'sorted');
    deg = degree(G);

    anchorCandidates = candidates(deg(candidates) >= 4);
    if isempty(anchorCandidates)
        anchorCandidates = candidates(deg(candidates) >= 3);
    end

    nAnchor = min(numel(anchorCandidates), max(0, round(anchorFraction * nToPick)));
    pickedAnchor = zeros(0, 1);
    if nAnchor > 0
        pickedAnchor = farthest_interior_sampling(G, anchorCandidates, seedIdx, nAnchor, nodeLabels);
    end

    remaining = setdiff(candidates, pickedAnchor, 'stable');
    nRemain = max(0, nToPick - numel(pickedAnchor));
    if nRemain <= 0 || isempty(remaining)
        picked = pickedAnchor;
        return;
    end

    nGrid = min(numel(remaining), round(0.45 * nRemain));
    pickedGrid = zeros(0, 1);
    if nGrid > 0
        pickedGrid = grid_interior_sampling( ...
            G, remaining, [seedIdx(:); pickedAnchor(:)], nGrid, nodeCoords, nodeLabels);
    end

    remaining2 = setdiff(remaining, pickedGrid, 'stable');
    nFar = max(0, nRemain - numel(pickedGrid));
    pickedFar = zeros(0, 1);
    if nFar > 0 && ~isempty(remaining2)
        pickedFar = farthest_interior_sampling( ...
            G, remaining2, [seedIdx(:); pickedAnchor(:); pickedGrid(:)], nFar, nodeLabels);
    end

    picked = unique([pickedAnchor(:); pickedGrid(:); pickedFar(:)], 'stable');
end

function picked = grid_interior_sampling(G, interiorIdx, seedIdx, nToPick, nodeCoords, nodeLabels)
    picked = zeros(0, 1);
    if nToPick <= 0 || isempty(interiorIdx)
        return;
    end

    candidates = unique(interiorIdx(:), 'sorted');
    nCand = numel(candidates);
    if nCand == 0
        return;
    end
    if nCand <= nToPick
        picked = candidates;
        return;
    end

    coords2d = nodeCoords(candidates, 1:min(2, size(nodeCoords, 2)));
    if size(coords2d, 2) == 1
        coords2d(:, 2) = 0;
    end
    x = coords2d(:, 1);
    y = coords2d(:, 2);

    xMin = min(x);
    xMax = max(x);
    yMin = min(y);
    yMax = max(y);
    xRange = xMax - xMin;
    yRange = yMax - yMin;

    if xRange < 1e-14 && yRange < 1e-14
        picked = farthest_interior_sampling(G, candidates, seedIdx, nToPick, nodeLabels);
        return;
    end

    nCells = max(1, min(nCand, round(nToPick)));
    if yRange < 1e-14
        nx = nCells;
        ny = 1;
    elseif xRange < 1e-14
        nx = 1;
        ny = nCells;
    else
        aspect = xRange / yRange;
        nx = max(1, round(sqrt(nCells * aspect)));
        ny = max(1, ceil(nCells / nx));
        while nx * ny < nCells
            nx = nx + 1;
        end
    end

    if xRange < 1e-14
        ix = ones(nCand, 1);
    else
        wx = xRange / nx;
        ix = floor((x - xMin) ./ wx) + 1;
        ix = min(max(ix, 1), nx);
    end

    if yRange < 1e-14
        iy = ones(nCand, 1);
    else
        wy = yRange / ny;
        iy = floor((y - yMin) ./ wy) + 1;
        iy = min(max(iy, 1), ny);
    end

    cellId = sub2ind([nx, ny], ix, iy);
    [uCell, ~, ic] = unique(cellId, 'sorted');
    reps = zeros(numel(uCell), 1);
    for c = 1:numel(uCell)
        members = find(ic == c);
        memberNodes = candidates(members);
        xy = coords2d(members, :);
        center = mean(xy, 1);
        d2 = sum((xy - center) .^ 2, 2);
        dMin = min(d2);
        tiePos = find(abs(d2 - dMin) <= 1e-15);
        if numel(tiePos) > 1
            [~, ord] = sort(nodeLabels(memberNodes(tiePos)), 'ascend');
            tiePos = tiePos(ord);
        end
        reps(c) = memberNodes(tiePos(1));
    end
    reps = unique(reps, 'stable');

    if numel(reps) > nToPick
        picked = farthest_interior_sampling(G, reps, seedIdx, nToPick, nodeLabels);
        return;
    end

    picked = reps;
    if numel(picked) < nToPick
        remaining = setdiff(candidates, picked, 'stable');
        nExtra = nToPick - numel(picked);
        if nExtra > 0 && ~isempty(remaining)
            extra = farthest_interior_sampling(G, remaining, [seedIdx(:); picked(:)], nExtra, nodeLabels);
            picked = unique([picked(:); extra(:)], 'stable');
        end
    end
end

function picked = farthest_interior_sampling(G, interiorIdx, seedIdx, nToPick, nodeLabels)
    picked = zeros(0, 1);
    if nToPick <= 0 || isempty(interiorIdx)
        return;
    end

    candidates = unique(interiorIdx(:), 'sorted');
    seedIdx = unique(seedIdx(:), 'sorted');

    if isempty(seedIdx)
        [~, k0] = min(nodeLabels(candidates));
        firstNode = candidates(k0);
        picked(end + 1, 1) = firstNode; %#ok<AGROW>
        candidates(k0) = [];
        seedIdx = firstNode;
    end

    if isempty(candidates) || numel(picked) >= nToPick
        return;
    end

    d0 = distances(G, candidates, seedIdx);
    if isvector(d0)
        minDist = d0(:);
    else
        minDist = min(d0, [], 2);
    end

    while ~isempty(candidates) && numel(picked) < nToPick
        maxDist = max(minDist);
        tiePos = find(minDist == maxDist);
        [~, order] = sort(nodeLabels(candidates(tiePos)), 'ascend');
        pickPos = tiePos(order(1));
        pickNode = candidates(pickPos);

        picked(end + 1, 1) = pickNode; %#ok<AGROW>
        candidates(pickPos) = [];
        minDist(pickPos) = [];

        if isempty(candidates)
            break;
        end

        dNew = distances(G, candidates, pickNode);
        minDist = min(minDist, dNew(:));
    end
end

function selectedOut = augment_connectivity_with_shortest_paths(G, selectedIn, nodeLabels)
    selectedOut = unique(selectedIn(:), 'sorted');
    if numel(selectedOut) <= 1
        return;
    end

    connectedSet = selectedOut(1);
    remaining = selectedOut(2:end);

    while ~isempty(remaining)
        sourceNode = remaining(1);
        d = distances(G, sourceNode, connectedSet);
        [dMin, idxMin] = min(d);

        if isinf(dMin)
            connectedSet = unique([connectedSet; sourceNode], 'sorted');
            remaining(1) = [];
            continue;
        end

        tiePos = find(d == dMin);
        if numel(tiePos) > 1
            targets = connectedSet(tiePos);
            [~, k] = min(nodeLabels(targets));
            targetNode = targets(k);
        else
            targetNode = connectedSet(idxMin);
        end

        pathNodes = shortestpath(G, sourceNode, targetNode);
        connectedSet = unique([connectedSet; pathNodes(:)], 'sorted');

        remaining = setdiff(remaining, connectedSet, 'stable');
    end

    selectedOut = connectedSet;
end

function edgesSubset = build_subset_edges(edgeListFull, subsetIdx, mapFullToSubset)
    edgesSubset = zeros(0, 2);
    if isempty(edgeListFull) || isempty(subsetIdx) || isempty(mapFullToSubset)
        return;
    end

    keep = mapFullToSubset(edgeListFull(:, 1)) > 0 & mapFullToSubset(edgeListFull(:, 2)) > 0;
    if ~any(keep)
        return;
    end

    e = edgeListFull(keep, :);
    edgesSubset = [mapFullToSubset(e(:, 1)), mapFullToSubset(e(:, 2))];
    edgesSubset = edgesSubset(edgesSubset(:, 1) > 0 & edgesSubset(:, 2) > 0, :);
    edgesSubset = edgesSubset(edgesSubset(:, 1) ~= edgesSubset(:, 2), :);
    if isempty(edgesSubset)
        edgesSubset = zeros(0, 2);
        return;
    end
    edgesSubset = sort(edgesSubset, 2);
    edgesSubset = unique(edgesSubset, 'rows', 'stable');
end

function boundaryEdges = get_boundary_edges(fullGraph, boundaryMask)
    boundaryEdges = zeros(0, 2);
    if isfield(fullGraph, 'boundary_edges_local') && ~isempty(fullGraph.boundary_edges_local)
        boundaryEdges = fullGraph.boundary_edges_local;
    elseif isfield(fullGraph, 'edges_local') && ~isempty(fullGraph.edges_local)
        e = fullGraph.edges_local;
        boundaryEdges = e(boundaryMask(e(:, 1)) & boundaryMask(e(:, 2)), :);
    end

    if ~isempty(boundaryEdges)
        boundaryEdges = sort(boundaryEdges, 2);
        boundaryEdges = unique(boundaryEdges, 'rows', 'stable');
    end
end

function boundaryEdgesSubset = map_boundary_edges_to_subset(fullGraph, subsetIdx, nFullNodes)
    boundaryEdgesSubset = zeros(0, 2);
    if isempty(subsetIdx)
        return;
    end

    boundaryEdgesFull = get_boundary_edges(fullGraph, fullGraph.boundary_mask(:));
    if isempty(boundaryEdgesFull)
        return;
    end

    mapFullToSubset = zeros(nFullNodes, 1);
    mapFullToSubset(subsetIdx) = 1:numel(subsetIdx);

    keep = mapFullToSubset(boundaryEdgesFull(:, 1)) > 0 & mapFullToSubset(boundaryEdgesFull(:, 2)) > 0;
    if ~any(keep)
        return;
    end

    b = boundaryEdgesFull(keep, :);
    boundaryEdgesSubset = [mapFullToSubset(b(:, 1)), mapFullToSubset(b(:, 2))];
    boundaryEdgesSubset = sort(boundaryEdgesSubset, 2);
    boundaryEdgesSubset = unique(boundaryEdgesSubset, 'rows', 'stable');
end

function keepRatio = derive_line_keep_ratio(detailLevel, lineNodeBudget, nReferenceNodes)
    if isempty(lineNodeBudget) || nReferenceNodes <= 0
        keepRatio = max(0.02, min(1, detailLevel));
        return;
    end
    keepRatio = min(1, max(0.01, lineNodeBudget / nReferenceNodes));
end

function [keptNodes, simplifiedEdges] = simplify_polyline_components(edgeList, nodeCoords, nodeLabels, keepRatio, tolDeg)
    keptNodes = zeros(0, 1);
    simplifiedEdges = zeros(0, 2);

    if isempty(edgeList)
        return;
    end

    nNodes = size(nodeCoords, 1);
    edgeList = sort(edgeList, 2);
    edgeList = unique(edgeList, 'rows', 'stable');
    edgeList = edgeList(edgeList(:, 1) > 0 & edgeList(:, 2) > 0, :);
    edgeList = edgeList(edgeList(:, 1) <= nNodes & edgeList(:, 2) <= nNodes, :);
    if isempty(edgeList)
        return;
    end

    G = graph(edgeList(:, 1), edgeList(:, 2), [], nNodes);
    deg = degree(G);
    activeNodes = find(deg > 0);
    if isempty(activeNodes)
        return;
    end

    comp = conncomp(G);
    compIds = unique(comp(activeNodes));
    compMinLabel = inf(size(compIds));
    for i = 1:numel(compIds)
        nodesComp = activeNodes(comp(activeNodes) == compIds(i));
        compMinLabel(i) = min(nodeLabels(nodesComp));
    end
    [~, compOrder] = sort(compMinLabel, 'ascend');
    compIds = compIds(compOrder);

    keptAll = zeros(0, 1);
    edgesAll = zeros(0, 2);

    for i = 1:numel(compIds)
        cid = compIds(i);
        nodesComp = activeNodes(comp(activeNodes) == cid);
        isCompEdge = ismember(edgeList(:, 1), nodesComp) & ismember(edgeList(:, 2), nodesComp);
        edgesComp = edgeList(isCompEdge, :);

        [orderedNodes, isCycle, ok] = order_component_as_path_or_cycle(edgesComp, nodeLabels);
        if ~ok
            keptComp = nodesComp(:);
            edgesCompOut = edgesComp;
        else
            keptComp = simplify_ordered_polyline( ...
                orderedNodes(:), isCycle, nodeCoords, nodeLabels, keepRatio, tolDeg);
            if numel(keptComp) <= 1
                edgesCompOut = zeros(0, 2);
            else
                edgesCompOut = [keptComp(1:end - 1), keptComp(2:end)];
                if isCycle && numel(keptComp) > 2
                    edgesCompOut(end + 1, :) = [keptComp(end), keptComp(1)]; %#ok<AGROW>
                end
            end
        end

        keptAll = [keptAll; keptComp(:)]; %#ok<AGROW>
        edgesAll = [edgesAll; edgesCompOut]; %#ok<AGROW>
    end

    keptNodes = unique(keptAll, 'sorted');
    if isempty(edgesAll)
        simplifiedEdges = zeros(0, 2);
        return;
    end

    keepNodeMask = false(nNodes, 1);
    keepNodeMask(keptNodes) = true;
    validEdges = keepNodeMask(edgesAll(:, 1)) & keepNodeMask(edgesAll(:, 2)) & ...
        (edgesAll(:, 1) ~= edgesAll(:, 2));
    edgesAll = edgesAll(validEdges, :);
    if isempty(edgesAll)
        simplifiedEdges = zeros(0, 2);
    else
        simplifiedEdges = sort(edgesAll, 2);
        simplifiedEdges = unique(simplifiedEdges, 'rows', 'stable');
    end
end

function [orderedNodes, isCycle, ok] = order_component_as_path_or_cycle(edgesComp, nodeLabels)
    orderedNodes = zeros(0, 1);
    isCycle = false;
    ok = false;

    if isempty(edgesComp)
        return;
    end

    nodes = unique(edgesComp(:), 'sorted');
    if numel(nodes) == 1
        orderedNodes = nodes;
        isCycle = false;
        ok = true;
        return;
    end

    map = zeros(max(nodes), 1);
    map(nodes) = 1:numel(nodes);
    localEdges = [map(edgesComp(:, 1)), map(edgesComp(:, 2))];
    Gc = graph(localEdges(:, 1), localEdges(:, 2), [], numel(nodes));
    d = degree(Gc);
    if any(d > 2)
        return;
    end

    labelsComp = nodeLabels(nodes);
    if all(d == 2)
        isCycle = true;
        [~, start] = min(labelsComp);
    else
        endpoints = find(d == 1);
        if numel(endpoints) ~= 2
            return;
        end
        [~, k] = min(labelsComp(endpoints));
        start = endpoints(k);
    end

    orderedLocal = zeros(numel(nodes), 1);
    visited = false(numel(nodes), 1);
    current = start;
    prev = 0;
    count = 0;

    while true
        count = count + 1;
        orderedLocal(count) = current;
        visited(current) = true;

        neigh = neighbors(Gc, current);
        if prev ~= 0
            neigh = neigh(neigh ~= prev);
        end

        if isempty(neigh)
            break;
        end

        if numel(neigh) > 1
            if prev == 0
                [~, kMin] = min(labelsComp(neigh));
                nextNode = neigh(kMin);
            else
                return;
            end
        else
            nextNode = neigh(1);
        end

        if visited(nextNode)
            if isCycle && nextNode == start && count == numel(nodes)
                break;
            else
                return;
            end
        end

        prev = current;
        current = nextNode;
    end

    if count ~= numel(nodes)
        return;
    end

    orderedNodes = nodes(orderedLocal(1:count));
    ok = true;
end

function kept = simplify_ordered_polyline(orderedNodes, isCycle, nodeCoords, nodeLabels, keepRatio, tolDeg)
    n = numel(orderedNodes);
    if n <= 2
        kept = orderedNodes(:);
        return;
    end

    keepRatio = max(0, min(1, keepRatio));
    tolRad = tolDeg * pi / 180;

    mandatory = false(n, 1);
    for i = 1:n
        if ~isCycle && (i == 1 || i == n)
            mandatory(i) = true;
            continue;
        end

        if i == 1
            iPrev = n;
            iNext = i + 1;
        elseif i == n
            iPrev = i - 1;
            iNext = 1;
        else
            iPrev = i - 1;
            iNext = i + 1;
        end

        pPrev = nodeCoords(orderedNodes(iPrev), :);
        pCurr = nodeCoords(orderedNodes(i), :);
        pNext = nodeCoords(orderedNodes(iNext), :);

        v1 = pPrev - pCurr;
        v2 = pNext - pCurr;
        n1 = norm(v1);
        n2 = norm(v2);
        if n1 < 1e-14 || n2 < 1e-14
            mandatory(i) = true;
            continue;
        end

        cosang = dot(v1, v2) / (n1 * n2);
        cosang = max(-1, min(1, cosang));
        angleVal = acos(cosang);
        bend = abs(pi - angleVal);
        mandatory(i) = bend > tolRad;
    end

    minKeep = 2;
    if isCycle
        minKeep = 3;
    end
    target = max(minKeep, ceil(keepRatio * n));
    target = min(n, target);

    keepIdx = find(mandatory);
    if numel(keepIdx) < minKeep
        [~, extraOrder] = sort(nodeLabels(orderedNodes), 'ascend');
        kNeed = minKeep - numel(keepIdx);
        addIdx = extraOrder(1:kNeed);
        keepIdx = unique([keepIdx; addIdx(:)], 'stable');
    end

    if numel(keepIdx) < target
        candidate = setdiff((1:n).', keepIdx, 'stable');
        need = target - numel(keepIdx);
        if ~isempty(candidate) && need > 0
            pos = unique(round(linspace(1, numel(candidate), need)));
            pos = max(1, min(numel(candidate), pos));
            keepIdx = [keepIdx; candidate(pos(:))];
        end
    end

    keepIdx = unique(keepIdx, 'sorted');
    if numel(keepIdx) > target
        pick = unique(round(linspace(1, numel(keepIdx), target)));
        keepIdx = keepIdx(pick);
    end

    if numel(keepIdx) < minKeep
        keepIdx = unique([keepIdx; 1; n], 'sorted');
    end

    kept = orderedNodes(keepIdx(:));
end

function v = get_default(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        v = s.(fieldName);
    else
        v = defaultValue;
    end
end
