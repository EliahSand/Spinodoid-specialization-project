function jlGraph = build_junction_ligament_graph(structGraph)
%BUILD_JUNCTION_LIGAMENT_GRAPH  Convert a skeleton-pixel graph into a
%   junction-ligament (JL) graph where nodes are junctions/endpoints and
%   edges are ligaments with geometric features.
%
%   Input:  structGraph with fields node_coords (Nx2), node_radius (Nx1),
%           edges_local (Ex2), graph (MATLAB graph object)
%   Output: jlGraph with node features, edge_index, edge_features, etc.

    G = structGraph.graph;
    n = structGraph.num_nodes;
    coords = structGraph.node_coords;
    radii  = structGraph.node_radius;
    deg = degree(G);

    if n == 0
        error('build_junction_ligament_graph:EmptyInput', 'Input graph has no nodes.');
    end

    % -------------------------------------------------------------------------
    % Identify key nodes: junctions (deg>2) and endpoints (deg==1)
    % degree-2 nodes are interior chain pixels and will be merged into edges.
    % -------------------------------------------------------------------------
    isKey = (deg > 2) | (deg == 1);
    keyNodes = find(isKey);
    if isempty(keyNodes)
        % Degenerate: every node has deg==2 (a pure cycle). Keep all.
        error('build_junction_ligament_graph:NoKeyNodes', ...
              'Graph has no junctions or endpoints (pure cycle). Not supported yet.');
    end

    % Edge lookup for fast O(1) neighbor traversal
    nEdge = size(structGraph.edges_local, 1);
    edgeLookup = sparse([structGraph.edges_local(:,1); structGraph.edges_local(:,2)], ...
                        [structGraph.edges_local(:,2); structGraph.edges_local(:,1)], ...
                        [(1:nEdge).'; (1:nEdge).'], n, n);
    visitedEdges = false(nEdge, 1);

    % -------------------------------------------------------------------------
    % Trace every chain between key nodes, each becomes a JL edge
    % -------------------------------------------------------------------------
    jlEdges = zeros(0, 2);
    jlEdgeFeatures = zeros(0, 7);
    localNodeForSkel = zeros(n, 1);  % local JL-node index for each skeleton node

    % Handle pure cycles (no key nodes) - but we already checked for that
    % Handle chains between key nodes
    for u = keyNodes(:).'
        nbrs = neighbors(G, u);
        for v = nbrs(:).'
            eIdx = edgeLookup(u, v);
            if eIdx == 0 || visitedEdges(eIdx)
                continue;
            end
            path = trace_chain(u, v, G, deg, edgeLookup, visitedEdges);
            pathEdgeIds = path_to_edge_ids(path, edgeLookup);
            visitedEdges(pathEdgeIds) = true;
            assign_edge(path);
        end
    end

    % Any remaining unvisited edges belong to isolated cycles
    freeEdges = find(~visitedEdges);
    while ~isempty(freeEdges)
        u = structGraph.edges_local(freeEdges(1), 1);
        v = structGraph.edges_local(freeEdges(1), 2);
        path = trace_chain(u, v, G, deg, edgeLookup, visitedEdges);
        pathEdgeIds = path_to_edge_ids(path, edgeLookup);
        visitedEdges(pathEdgeIds) = true;
        assign_edge(path);
        freeEdges = find(~visitedEdges);
    end

    % -------------------------------------------------------------------------
    % Build JL node coordinates and features
    % -------------------------------------------------------------------------
    jlCoords = coords(keyNodes, :);
    jlRadii  = radii(keyNodes);
    jlDeg    = deg(keyNodes);
    nJL = numel(keyNodes);

    % Build edge_index in [2 x E] format
    E = size(jlEdges, 1);
    edge_index = jlEdges.';
    if isempty(edge_index), edge_index = zeros(2, 0, 'int32'); end

    jlGraph = struct();
    jlGraph.num_nodes = nJL;
    jlGraph.num_edges = E;
    jlGraph.node_coords = jlCoords;
    jlGraph.node_radius = jlRadii;
    jlGraph.node_degree = jlDeg;
    jlGraph.edges_local = jlEdges;
    jlGraph.edge_index = edge_index;
    jlGraph.edge_features = jlEdgeFeatures;
    jlGraph.skel_to_jl = keyNodes;  % mapping: keyNodes(i) is original skeleton node i

    % Build compact graph object for convenience
    if E == 0
        jlGraph.graph = graph([], [], [], nJL);
    else
        jlGraph.graph = graph(jlEdges(:,1), jlEdges(:,2), jlEdgeFeatures(:,1), nJL);
    end

    %% -----------------------------------------------------------------------
    % Local helpers
    %% -----------------------------------------------------------------------

    function assign_edge(path)
        if numel(path) < 2
            return;
        end
        srcSkel = path(1);
        dstSkel = path(end);
        % Map to JL node indices
        src = find(keyNodes == srcSkel, 1);
        dst = find(keyNodes == dstSkel, 1);
        if isempty(src) || isempty(dst)
            return;  % should not happen
        end
        if src == dst
            return;  % skip self-loops
        end

        chainCoords = coords(path, :);
        chainRadii  = radii(path);

        % Geometric features
        physLength = path_length(chainCoords);
        eucDist = norm(chainCoords(1,:) - chainCoords(end,:));
        avgR = mean(chainRadii, 'omitnan');
        stdR = std(chainRadii, 'omitnan');
        if ~isfinite(stdR), stdR = 0; end
        cvR = safe_div(stdR, avgR);
        tortuosity = safe_div(physLength, max(eucDist, eps));
        maxCurv = max_chain_curvature(chainCoords);
        orientVec = chainCoords(end,:) - chainCoords(1,:);
        orientVec = orientVec / max(norm(orientVec), eps);

        jlEdges(end+1, :) = [src, dst];
        jlEdgeFeatures(end+1, :) = [physLength, avgR, cvR, maxCurv, tortuosity, orientVec(1), orientVec(2)];
    end

end

%% -------------------------------------------------------------------------
% Trace a chain from startNode through nextNode until hitting a key node
%% -------------------------------------------------------------------------
function path = trace_chain(startNode, nextNode, G, deg, edgeLookup, visited)
    path = [startNode, nextNode];
    prev = startNode;
    curr = nextNode;
    while true
        if curr == startNode && numel(path) > 2
            break;  % closed cycle
        end
        if deg(curr) ~= 2
            break;  % hit a key node
        end
        nbrs = sort(neighbors(G, curr)).';
        nbrs(nbrs == prev) = [];
        if isempty(nbrs)
            break;
        end
        nxt = nbrs(1);
        eIdx = edgeLookup(curr, nxt);
        if eIdx == 0
            break;
        end
        if visited(eIdx) && nxt ~= startNode
            break;
        end
        path(end+1) = nxt;
        prev = curr;
        curr = nxt;
    end
end

function edgeIds = path_to_edge_ids(path, edgeLookup)
    edgeIds = zeros(numel(path)-1, 1);
    for i = 1:(numel(path)-1)
        e = edgeLookup(path(i), path(i+1));
        edgeIds(i) = e;
    end
    edgeIds = edgeIds(edgeIds > 0);
end

function len = path_length(coords)
    if size(coords, 1) < 2
        len = 0;
        return;
    end
    d = diff(coords, 1, 1);
    len = sum(sqrt(sum(d.^2, 2)));
end

function c = max_chain_curvature(coords)
    N = size(coords, 1);
    if N < 3
        c = 0;
        return;
    end
    curvs = zeros(N - 2, 1);
    for i = 2:(N - 1)
        a = coords(i-1, :) - coords(i, :);
        b = coords(i+1, :) - coords(i, :);
        an = norm(a);
        bn = norm(b);
        if an < eps || bn < eps
            curvs(i-1) = 0;
            continue;
        end
        cosTheta = dot(a, b) / (an * bn);
        cosTheta = max(-1, min(1, cosTheta));
        curvs(i-1) = acos(cosTheta);
    end
    c = max(curvs);
end

function out = safe_div(a, b)
    if abs(b) < eps
        out = 0;
    else
        out = a / b;
    end
end
