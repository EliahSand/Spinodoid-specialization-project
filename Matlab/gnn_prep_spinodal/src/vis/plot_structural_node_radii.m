function ax = plot_structural_node_radii(fullGraph, structGraph, varargin)
%PLOT_STRUCTURAL_NODE_RADII Draw structural nodes as radius-scaled disks.
%
% ax = PLOT_STRUCTURAL_NODE_RADII(fullGraph, structGraph, ...)
%
% The dense reference graph is drawn as a readable gray spinodoid
% background. Each structural node radius is drawn as a faint blue outline.
%
% Name-Value options:
%   'AxesHandle'  : target axes (default create new figure/axes)
%   'Title'       : plot title (default 'Structural Node Radii')
%   'UnitsScale'  : coordinate / radius scale factor (default 1)
%   'XLabel'      : x-axis label
%   'YLabel'      : y-axis label
%   'ShowDense'   : draw dense reference edges (default true)
%   'ShowEdges'   : draw structural graph edges (default true)

    p = inputParser;
    p.addParameter('AxesHandle', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
    p.addParameter('Title', 'Structural Node Radii', @(x) ischar(x) || isstring(x));
    p.addParameter('UnitsScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('XLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('YLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ShowDense', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('ShowEdges', true, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    denseEdgeColor = [0.52, 0.52, 0.52];
    denseEdgeWidth = 0.45;
    structuralEdgeColor = [0.45, 0.45, 0.45];
    structuralEdgeWidth = 0.45;

    if isempty(opts.AxesHandle)
        figure('Name', 'Structural Node Radii');
        ax = axes;
    else
        ax = opts.AxesHandle;
    end

    hold(ax, 'on');

    if logical(opts.ShowDense) && isfield(fullGraph, 'edges_local') && ~isempty(fullGraph.edges_local)
        draw_edges_local(ax, fullGraph.node_coords, fullGraph.edges_local, ...
            opts.UnitsScale, denseEdgeColor, denseEdgeWidth);
    end

    if logical(opts.ShowEdges) && isfield(structGraph, 'edges_local') && ~isempty(structGraph.edges_local)
        draw_edges_local(ax, structGraph.node_coords, structGraph.edges_local, ...
            opts.UnitsScale, structuralEdgeColor, structuralEdgeWidth);
    end

    xy = structGraph.node_coords(:, 1:2) * opts.UnitsScale;
    r = get_node_radii(structGraph);

    if ~isempty(r)
        r = r * opts.UnitsScale;
        valid = isfinite(r) & r > 0;
        if any(valid)
            draw_radius_disks(ax, xy(valid, 1), xy(valid, 2), r(valid));
        end
        scatter(ax, xy(:, 1), xy(:, 2), 10, [0.18, 0.18, 0.18], 'filled', 'MarkerEdgeColor', 'none');
    else
        scatter(ax, xy(:, 1), xy(:, 2), 20, [0.85, 0.20, 0.10], 'filled', 'MarkerEdgeColor', 'none');
    end

    axis(ax, 'equal');
    axis(ax, 'tight');
    grid(ax, 'off');
    box(ax, 'on');
    set(ax, 'FontSize', 14, 'LineWidth', 1.0);
    if strlength(string(opts.Title)) > 0
        title(ax, char(string(opts.Title)));
    end
    if strlength(string(opts.XLabel)) > 0
        xlabel(ax, char(string(opts.XLabel)));
    end
    if strlength(string(opts.YLabel)) > 0
        ylabel(ax, char(string(opts.YLabel)));
    end

    hold(ax, 'off');
end

function draw_edges_local(ax, coords, edges, unitsScale, color, lw)
    xy = coords(:, 1:2) * unitsScale;
    x = [xy(edges(:,1), 1), xy(edges(:,2), 1), nan(size(edges,1),1)].';
    y = [xy(edges(:,1), 2), xy(edges(:,2), 2), nan(size(edges,1),1)].';
    plot(ax, x(:), y(:), '-', 'Color', color, 'LineWidth', lw);
end

function r = get_node_radii(structGraph)
    r = [];
    if isfield(structGraph, 'node_radius') && numel(structGraph.node_radius) == size(structGraph.node_coords, 1)
        r = double(structGraph.node_radius(:));
        return;
    end
    if isfield(structGraph, 'node_data_table') && istable(structGraph.node_data_table) ...
            && any(strcmp(structGraph.node_data_table.Properties.VariableNames, 'Radius'))
        r = double(structGraph.node_data_table.Radius(:));
    end
end
