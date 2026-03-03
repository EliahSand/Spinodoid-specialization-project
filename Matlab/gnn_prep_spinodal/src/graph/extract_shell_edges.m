function [edges, ownerElementIdx] = extract_shell_edges(elementLocalNodes)
%EXTRACT_SHELL_EDGES Convert element connectivity to undirected shell edges.
%
% [edges, ownerElementIdx] = EXTRACT_SHELL_EDGES(elementLocalNodes)
%
% Input:
%   elementLocalNodes - Mx1 cell; each cell is ordered node indices of one element.
%
% Output:
%   edges           - Px2 undirected node-index pairs (not unique, includes duplicates)
%   ownerElementIdx - Px1 element index for each edge row in edges

    nElem = numel(elementLocalNodes);
    edgeBlocks = cell(nElem, 1);
    ownerBlocks = cell(nElem, 1);

    for e = 1:nElem
        nodes = round(elementLocalNodes{e}(:));
        nodes = nodes(nodes > 0);
        if numel(nodes) < 2
            edgeBlocks{e} = zeros(0, 2);
            ownerBlocks{e} = zeros(0, 1);
            continue;
        end

        nextNodes = nodes([2:end, 1]);
        ePairs = [nodes, nextNodes];
        ePairs = sort(ePairs, 2);
        ePairs = ePairs(ePairs(:, 1) ~= ePairs(:, 2), :);

        edgeBlocks{e} = ePairs;
        ownerBlocks{e} = repmat(e, size(ePairs, 1), 1);
    end

    if isempty(edgeBlocks)
        edges = zeros(0, 2);
        ownerElementIdx = zeros(0, 1);
        return;
    end

    edges = vertcat(edgeBlocks{:});
    ownerElementIdx = vertcat(ownerBlocks{:});
end

