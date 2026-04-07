function ax = plot_structural_graph_cleanup_overlay(fullGraph, skeletonGraph, structGraphBefore, structGraphAfter, cleanupDebug, varargin)
%PLOT_STRUCTURAL_GRAPH_CLEANUP_OVERLAY Compare structural graph before/after cleanup.
%
% ax = PLOT_STRUCTURAL_GRAPH_CLEANUP_OVERLAY(fullGraph, skeletonGraph, structGraphBefore, structGraphAfter, cleanupDebug, ...)
%
% Name-Value options:
%   'AxesHandle'   : target axes (default create new figure/axes)
%   'Title'        : plot title
%   'ShowDense'    : show dense reference graph (default true)
%   'ShowSkeleton' : show skeleton graph (default true)
%   'UnitsScale'   : coordinate scale factor (default 1)
%   'XLabel'       : x-axis label
%   'YLabel'       : y-axis label

    p = inputParser;
    p.addParameter('AxesHandle', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
    p.addParameter('Title', 'Structural Graph Cleanup Overlay', @(x) ischar(x) || isstring(x));
    p.addParameter('ShowDense', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ShowSkeleton', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('UnitsScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('XLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('YLabel', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    if isempty(opts.AxesHandle)
        fig = figure('Name', 'Structural Graph Cleanup Overlay'); %#ok<NASGU>
        ax = axes;
    else
        ax = opts.AxesHandle;
    end

    hold(ax, 'on');
    legendHandles = gobjects(0);
    legendLabels = {};

    if logical(opts.ShowDense)
        h = draw_graph_edges(ax, fullGraph.node_coords, fullGraph.edges_local, opts.UnitsScale, [0.82, 0.82, 0.82], 0.30, '-');
        if ~isempty(h)
            legendHandles(end + 1) = h; %#ok<AGROW>
            legendLabels{end + 1} = 'Dense'; %#ok<AGROW>
        end
    end

    if logical(opts.ShowSkeleton)
        h = draw_graph_edges(ax, skeletonGraph.node_coords, skeletonGraph.edges_local, opts.UnitsScale, [0.20, 0.52, 0.80], 0.90, '-');
        if ~isempty(h)
            legendHandles(end + 1) = h; %#ok<AGROW>
            legendLabels{end + 1} = 'Skeleton'; %#ok<AGROW>
        end
    end

    h = draw_graph_edges(ax, structGraphBefore.node_coords, structGraphBefore.edges_local, opts.UnitsScale, [0.55, 0.55, 0.55], 1.10, '--');
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural before cleanup'; %#ok<AGROW>
    end
    h = draw_graph_nodes(ax, structGraphBefore.node_coords, opts.UnitsScale, [0.45, 0.45, 0.45], 18, 'o');
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Before nodes'; %#ok<AGROW>
    end

    h = draw_graph_edges(ax, structGraphAfter.node_coords, structGraphAfter.edges_local, opts.UnitsScale, [0.92, 0.42, 0.12], 1.8, '-');
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural after cleanup'; %#ok<AGROW>
    end
    h = draw_graph_nodes(ax, structGraphAfter.node_coords, opts.UnitsScale, [0.92, 0.22, 0.12], 26, 'o');
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'After nodes'; %#ok<AGROW>
    end

    removedCoords = zeros(0, 2);
    if nargin >= 5 && isstruct(cleanupDebug) && isfield(cleanupDebug, 'removed_node_coords')
        removedCoords = cleanupDebug.removed_node_coords;
    end
    if ~isempty(removedCoords)
        h = draw_graph_nodes(ax, removedCoords, opts.UnitsScale, [0.78, 0.12, 0.65], 40, 'x');
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Removed nodes'; %#ok<AGROW>
    end

    axis(ax, 'equal');
    box(ax, 'on');
    grid(ax, 'on');
    if strlength(string(opts.Title)) > 0
        title(ax, char(string(opts.Title)));
    end
    if strlength(string(opts.XLabel)) > 0
        xlabel(ax, char(string(opts.XLabel)));
    end
    if strlength(string(opts.YLabel)) > 0
        ylabel(ax, char(string(opts.YLabel)));
    end

    if ~isempty(legendHandles)
        legend(ax, legendHandles, legendLabels, 'Location', 'bestoutside');
    end
end

function h = draw_graph_edges(ax, coords, edges, unitsScale, color, lineWidth, lineStyle)
    if isempty(edges)
        h = gobjects(0);
        return;
    end
    xy = ensure_xy(coords) * unitsScale;
    x = [xy(edges(:, 1), 1), xy(edges(:, 2), 1), nan(size(edges, 1), 1)].';
    y = [xy(edges(:, 1), 2), xy(edges(:, 2), 2), nan(size(edges, 1), 1)].';
    h = plot(ax, x(:), y(:), 'LineStyle', lineStyle, 'Color', color, 'LineWidth', lineWidth);
end

function h = draw_graph_nodes(ax, coords, unitsScale, color, markerSize, marker)
    if isempty(coords)
        h = gobjects(0);
        return;
    end
    xy = ensure_xy(coords) * unitsScale;
    h = scatter(ax, xy(:, 1), xy(:, 2), markerSize, ...
        'Marker', marker, ...
        'MarkerEdgeColor', color, ...
        'LineWidth', 1.1);
    if ~strcmp(marker, 'x')
        set(h, 'MarkerFaceColor', color);
    end
end

function xy = ensure_xy(coords)
    if size(coords, 2) < 2
        error('plot_structural_graph_cleanup_overlay:BadCoords', ...
            'node_coords must have at least two columns [X,Y].');
    end
    xy = coords(:, 1:2);
end
