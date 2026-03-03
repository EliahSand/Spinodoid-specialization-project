function fullGraph = build_full_reference_graph(inpData, varargin)
%BUILD_FULL_REFERENCE_GRAPH Build full graph for selected top/spinodal region.
%
% fullGraph = BUILD_FULL_REFERENCE_GRAPH(inpData, ...)
%
% Name-Value options:
%   'ElsetName'         : explicit elset name (default '')
%   'AutoDetectElset'   : auto-detect spinodal/top elset when ElsetName empty
%   'PatternPriority'   : patterns used by auto-detect
%
% Output fields:
%   node_labels         - graph node labels (Abaqus labels)
%   node_coords         - graph node coordinates
%   edges_local         - Ex2 unique undirected edges (local node indices)
%   boundary_mask       - logical mask over graph nodes
%   boundary_labels     - Abaqus labels of boundary nodes
%   element_labels      - selected region element labels
%   element_node_labels - selected region element connectivity (Abaqus labels)
%   element_local_nodes - selected region element connectivity (local indices)
%   graph               - MATLAB graph object

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
        error('build_full_reference_graph:EmptyRegion', ...
            'Selected region has no nodes. Check elset selection.');
    end

    [isNodeKnown, nodeCoordIdx] = ismember(nodeLabels, inpData.node_labels);
    if any(~isNodeKnown)
        error('build_full_reference_graph:MissingNodeCoords', ...
            'Some selected nodes are missing from node coordinate table.');
    end
    nodeCoords = inpData.node_coords(nodeCoordIdx, :);

    nElem = numel(region.element_connectivity);
    elementLocalNodes = cell(nElem, 1);
    for e = 1:nElem
        elemNodeLabels = region.element_connectivity{e}(:);
        [ok, localIdx] = ismember(elemNodeLabels, nodeLabels);
        if any(~ok)
            error('build_full_reference_graph:ElementMappingFailed', ...
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
    fullGraph.kind = 'full_reference_graph';
    fullGraph.source_inp_path = inpData.path;
    fullGraph.region_source = region.source;
    fullGraph.region_elset_name = region.elset_name;

    fullGraph.node_labels = nodeLabels;
    fullGraph.node_coords = nodeCoords;
    fullGraph.num_nodes = numel(nodeLabels);

    fullGraph.element_labels = region.element_labels(:);
    fullGraph.element_node_labels = region.element_connectivity;
    fullGraph.element_local_nodes = elementLocalNodes;
    fullGraph.num_elements = nElem;

    fullGraph.edges_local = edgesLocal;
    fullGraph.edge_use_count = boundary.edge_use_count;
    fullGraph.boundary_edges_local = boundary.boundary_edges;
    fullGraph.boundary_mask = boundary.boundary_mask;
    fullGraph.boundary_labels = nodeLabels(boundary.boundary_mask);

    fullGraph.graph = G;
    fullGraph.node_index_table = table((1:fullGraph.num_nodes).', nodeLabels, ...
        'VariableNames', {'NodeIndex', 'Label'});
end

