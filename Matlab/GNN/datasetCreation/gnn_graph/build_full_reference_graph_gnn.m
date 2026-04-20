function fullGraph = build_full_reference_graph_gnn(inpData, varargin)
%BUILD_FULL_REFERENCE_GRAPH_GNN Minimal variant of build_full_reference_graph for GNN export.
%
% Produces only the fields required by extract_structural_graph_gnn:
%   node_labels   - Abaqus node labels
%   node_coords   - (N x 3) node coordinates
%   num_nodes     - scalar
%   edges_local   - (E x 2) unique undirected edge pairs (local indices)
%   boundary_mask - (N x 1) logical, true for boundary nodes
%   graph         - MATLAB graph object
%
% All element connectivity, region metadata, and index tables are omitted.

    p = inputParser;
    p.addParameter('ElsetName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoDetectElset', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PatternPriority', {'spinodal', 'top'}, @(x) iscell(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    region = select_region_elements(inpData, ...
        'ElsetName', opts.ElsetName, ...
        'AutoDetectElset', logical(opts.AutoDetectElset), ...
        'PatternPriority', opts.PatternPriority);

    nodeLabels = region.node_labels(:);
    if isempty(nodeLabels)
        error('build_full_reference_graph_gnn:EmptyRegion', ...
            'Selected region has no nodes. Check elset selection.');
    end

    [isNodeKnown, nodeCoordIdx] = ismember(nodeLabels, inpData.node_labels);
    if any(~isNodeKnown)
        error('build_full_reference_graph_gnn:MissingNodeCoords', ...
            'Some selected nodes are missing from node coordinate table.');
    end
    nodeCoords = inpData.node_coords(nodeCoordIdx, :);

    nElem = numel(region.element_connectivity);
    elementLocalNodes = cell(nElem, 1);
    for e = 1:nElem
        elemNodeLabels = region.element_connectivity{e}(:);
        [ok, localIdx] = ismember(elemNodeLabels, nodeLabels);
        if any(~ok)
            error('build_full_reference_graph_gnn:ElementMappingFailed', ...
                'Element %d references node labels outside selected node set.', ...
                region.element_labels(e));
        end
        elementLocalNodes{e} = localIdx(:);
    end

    boundary = detect_boundary_nodes(elementLocalNodes, numel(nodeLabels));
    edgesLocal = boundary.unique_edges;
    if isempty(edgesLocal)
        G = graph([], [], [], numel(nodeLabels));
    else
        G = graph(edgesLocal(:, 1), edgesLocal(:, 2), [], numel(nodeLabels));
    end

    fullGraph = struct();
    fullGraph.node_labels = nodeLabels;
    fullGraph.node_coords = nodeCoords;
    fullGraph.num_nodes = numel(nodeLabels);
    fullGraph.edges_local = edgesLocal;
    fullGraph.boundary_mask = boundary.boundary_mask;
    fullGraph.graph = G;
end
