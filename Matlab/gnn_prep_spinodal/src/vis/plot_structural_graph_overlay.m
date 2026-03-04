function ax = plot_structural_graph_overlay(fullGraph, skeletonGraph, structGraph, varargin)
%PLOT_STRUCTURAL_GRAPH_OVERLAY Overlay dense, skeleton, and reduced structural graphs.
%
% ax = PLOT_STRUCTURAL_GRAPH_OVERLAY(fullGraph, skeletonGraph, structGraph, ...)
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
    p.addParameter('Title', 'Structural Graph Overlay', @(x) ischar(x) || isstring(x));
    p.addParameter('ShowDense', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ShowSkeleton', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('UnitsScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('XLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('YLabel', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    if isempty(opts.AxesHandle)
        fig = figure('Name', 'Structural Graph Overlay'); %#ok<NASGU>
        ax = axes;
    else
        ax = opts.AxesHandle;
    end

    hold(ax, 'on');

    legendHandles = gobjects(0);
    legendLabels = {};

    if logical(opts.ShowDense)
        h = draw_graph_edges(ax, fullGraph.node_coords, fullGraph.edges_local, opts.UnitsScale, [0.75, 0.75, 0.75], 0.35);
        if ~isempty(h)
            legendHandles(end + 1) = h; %#ok<AGROW>
            legendLabels{end + 1} = 'Dense'; %#ok<AGROW>
        end
    end

    if logical(opts.ShowSkeleton)
        h = draw_graph_edges(ax, skeletonGraph.node_coords, skeletonGraph.edges_local, opts.UnitsScale, [0.18, 0.55, 0.82], 1.05);
        if ~isempty(h)
            legendHandles(end + 1) = h; %#ok<AGROW>
            legendLabels{end + 1} = 'Skeleton'; %#ok<AGROW>
        end
        h = draw_graph_nodes(ax, skeletonGraph.node_coords, opts.UnitsScale, [0.18, 0.55, 0.82], 10);
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Skeleton nodes'; %#ok<AGROW>
    end

    h = draw_graph_edges(ax, structGraph.node_coords, structGraph.edges_local, opts.UnitsScale, [0.92, 0.42, 0.12], 1.8);
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural'; %#ok<AGROW>
    end
    h = draw_graph_nodes(ax, structGraph.node_coords, opts.UnitsScale, [0.92, 0.22, 0.12], 24);
    legendHandles(end + 1) = h; %#ok<AGROW>
    legendLabels{end + 1} = 'Structural nodes'; %#ok<AGROW>

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

function h = draw_graph_edges(ax, coords, edges, unitsScale, color, lineWidth)
    if isempty(edges)
        h = gobjects(0);
        return;
    end
    xy = ensure_xy(coords) * unitsScale;
    x = [xy(edges(:, 1), 1), xy(edges(:, 2), 1), nan(size(edges, 1), 1)].';
    y = [xy(edges(:, 1), 2), xy(edges(:, 2), 2), nan(size(edges, 1), 1)].';
    h = plot(ax, x(:), y(:), '-', 'Color', color, 'LineWidth', lineWidth);
end

function h = draw_graph_nodes(ax, coords, unitsScale, color, markerSize)
    xy = ensure_xy(coords) * unitsScale;
    h = scatter(ax, xy(:, 1), xy(:, 2), markerSize, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
end

function xy = ensure_xy(coords)
    if size(coords, 2) < 2
        error('plot_structural_graph_overlay:BadCoords', ...
            'node_coords must have at least two columns [X,Y].');
    end
    xy = coords(:, 1:2);
end
