function ax = plot_structural_thickness_overlay(fullGraph, structGraph, varargin)
%PLOT_STRUCTURAL_THICKNESS_OVERLAY Overlay dense graph with structural-node thickness coloring.
%
% ax = PLOT_STRUCTURAL_THICKNESS_OVERLAY(fullGraph, structGraph, ...)
%
% Name-Value options:
%   'AxesHandle' : target axes (default create new figure/axes)
%   'Title'      : plot title
%   'UnitsScale' : coordinate and thickness scale factor (default 1)
%   'XLabel'     : x-axis label
%   'YLabel'     : y-axis label

    p = inputParser;
    p.addParameter('AxesHandle', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
    p.addParameter('Title', 'Structural Thickness Overlay', @(x) ischar(x) || isstring(x));
    p.addParameter('UnitsScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('XLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('YLabel', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    if isempty(opts.AxesHandle)
        fig = figure('Name', 'Structural Thickness Overlay'); %#ok<NASGU>
        ax = axes;
    else
        ax = opts.AxesHandle;
    end

    hold(ax, 'on');
    colormap(ax, parula);

    legendHandles = gobjects(0);
    legendLabels = {};

    h = draw_graph_edges(ax, fullGraph.node_coords, fullGraph.edges_local, opts.UnitsScale, [0.82, 0.82, 0.82], 0.30);
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Dense'; %#ok<AGROW>
    end

    h = draw_graph_edges(ax, structGraph.node_coords, structGraph.edges_local, opts.UnitsScale, [0.25, 0.25, 0.25], 1.30);
    if ~isempty(h)
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural'; %#ok<AGROW>
    end

    thicknessVals = structural_node_thickness(structGraph);
    if isempty(thicknessVals)
        h = draw_graph_nodes_plain(ax, structGraph.node_coords, opts.UnitsScale, [0.92, 0.22, 0.12], 28);
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural nodes'; %#ok<AGROW>
    else
        h = draw_graph_nodes_colored(ax, structGraph.node_coords, thicknessVals * opts.UnitsScale, opts.UnitsScale, 30);
        legendHandles(end + 1) = h; %#ok<AGROW>
        legendLabels{end + 1} = 'Structural nodes (thickness)'; %#ok<AGROW>
        cb = colorbar(ax);
        cb.Label.String = thickness_label(opts.UnitsScale);
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

function h = draw_graph_nodes_plain(ax, coords, unitsScale, color, markerSize)
    xy = ensure_xy(coords) * unitsScale;
    h = scatter(ax, xy(:, 1), xy(:, 2), markerSize, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
end

function h = draw_graph_nodes_colored(ax, coords, values, unitsScale, markerSize)
    xy = ensure_xy(coords) * unitsScale;
    h = scatter(ax, xy(:, 1), xy(:, 2), markerSize, values, 'filled', ...
        'MarkerEdgeColor', [0.20, 0.20, 0.20], 'LineWidth', 0.25);
end

function vals = structural_node_thickness(structGraph)
    vals = [];
    if isfield(structGraph, 'node_thickness') && numel(structGraph.node_thickness) == size(structGraph.node_coords, 1)
        vals = double(structGraph.node_thickness(:));
        return;
    end
    if isfield(structGraph, 'node_data_table') && istable(structGraph.node_data_table) ...
            && any(strcmp(structGraph.node_data_table.Properties.VariableNames, 'Thickness'))
        vals = double(structGraph.node_data_table.Thickness(:));
    end
end

function txt = thickness_label(unitsScale)
    if abs(unitsScale - 1000) < 1e-12
        txt = 'Thickness [mm]';
    elseif abs(unitsScale - 1) < 1e-12
        txt = 'Thickness [m]';
    else
        txt = 'Thickness [scaled]';
    end
end

function xy = ensure_xy(coords)
    if size(coords, 2) < 2
        error('plot_structural_thickness_overlay:BadCoords', ...
            'node_coords must have at least two columns [X,Y].');
    end
    xy = coords(:, 1:2);
end
