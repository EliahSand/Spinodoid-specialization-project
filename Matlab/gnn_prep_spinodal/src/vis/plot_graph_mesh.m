function h = plot_graph_mesh(graphData, varargin)
%PLOT_GRAPH_MESH Plot graph nodes/edges with optional boundary highlight.
%
% h = PLOT_GRAPH_MESH(graphData, ...)
%
% Name-Value options:
%   'Title'         : plot title (default '')
%   'AxesHandle'    : target axes (default create new figure/axes)
%   'ShowBoundary'  : legacy boolean boundary toggle (default true)
%   'BoundaryMode'  : 'all' | 'interior' | 'none' (default 'all')
%   'Style'         : 'default' | 'paper2d' (default 'default')
%   'UnitsScale'    : multiply coordinates before plotting (default 1)
%   'XLabel'        : x-axis label (default '')
%   'YLabel'        : y-axis label (default '')
%   'NodeMarkerSize'     : marker size for regular nodes (optional)
%   'BoundaryMarkerSize' : marker size for highlighted boundary nodes (optional)

    p = inputParser;
    p.addParameter('Title', '', @(x) ischar(x) || isstring(x));
    p.addParameter('AxesHandle', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
    p.addParameter('ShowBoundary', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('BoundaryMode', 'all', @(x) ischar(x) || isstring(x));
    p.addParameter('Style', 'default', @(x) ischar(x) || isstring(x));
    p.addParameter('UnitsScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('XLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('YLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('NodeMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.addParameter('BoundaryMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.parse(varargin{:});
    opts = p.Results;
    style = lower(char(string(opts.Style)));
    style = validatestring(style, {'default', 'paper2d'});
    boundaryMode = lower(char(string(opts.BoundaryMode)));
    boundaryMode = validatestring(boundaryMode, {'all', 'interior', 'none'});
    if ~logical(opts.ShowBoundary)
        boundaryMode = 'none';
    end

    if isempty(opts.AxesHandle)
        fig = figure;
        ax = axes;
    else
        fig = ancestor(opts.AxesHandle, 'figure');
        ax = opts.AxesHandle;
    end

    coords = graphData.node_coords;
    if size(coords, 2) < 2
        error('plot_graph_mesh:BadCoords', ...
            'graphData.node_coords must have at least two columns [X,Y].');
    end
    if size(coords, 2) < 3
        coords(:, 3) = 0;
    end
    x = coords(:, 1) * opts.UnitsScale;
    y = coords(:, 2) * opts.UnitsScale;
    z = coords(:, 3) * opts.UnitsScale;

    if isfield(graphData, 'graph') && isa(graphData.graph, 'graph')
        G = graphData.graph;
    else
        edges = graphData.edges_local;
        G = graph(edges(:, 1), edges(:, 2), [], size(coords, 1));
    end

    is2D = all(abs(z - z(1)) < 1e-14);
    if is2D
        h = plot(ax, G, 'XData', x, 'YData', y, 'NodeLabel', {});
    else
        h = plot(ax, G, 'XData', x, 'YData', y, 'ZData', z, 'NodeLabel', {});
        view(ax, 3);
    end

    switch style
        case 'paper2d'
            h.NodeColor = [0.18, 0.42, 0.85];
            h.EdgeColor = [0.18, 0.42, 0.85];
            h.LineWidth = 0.9;
            h.Marker = 'o';
            if isempty(opts.NodeMarkerSize)
                h.MarkerSize = 3;
            else
                h.MarkerSize = opts.NodeMarkerSize;
            end
            ax.Color = [0.93, 0.93, 0.93];
            if ~isempty(fig)
                fig.Color = [0.93, 0.93, 0.93];
            end
            grid(ax, 'off');
            box(ax, 'on');
            set(ax, 'FontSize', 18, 'LineWidth', 1.0);

        otherwise
            h.NodeColor = [0.10, 0.25, 0.65];
            h.EdgeColor = [0.35, 0.35, 0.35];
            h.Marker = '.';
            if isempty(opts.NodeMarkerSize)
                h.MarkerSize = 8;
            else
                h.MarkerSize = opts.NodeMarkerSize;
            end
            grid(ax, 'on');
    end

    boundaryMarkerSize = opts.BoundaryMarkerSize;
    if isempty(boundaryMarkerSize)
        if strcmp(style, 'paper2d')
            boundaryMarkerSize = 4.5;
        else
            boundaryMarkerSize = max(h.MarkerSize + 1, 6);
        end
    end

    if ~strcmp(boundaryMode, 'none') && isfield(graphData, 'boundary_mask')
        bNodes = find(logical(graphData.boundary_mask(:)));
        if strcmp(boundaryMode, 'interior')
            bNodes = interior_boundary_nodes([x, y, z], bNodes);
        end
        if ~isempty(bNodes)
            highlight(h, bNodes, 'NodeColor', [0.85, 0.20, 0.20], ...
                'MarkerSize', boundaryMarkerSize);
        end
    end

    axis(ax, 'equal');
    if strlength(string(opts.XLabel)) > 0
        xlabel(ax, char(string(opts.XLabel)));
    end
    if strlength(string(opts.YLabel)) > 0
        ylabel(ax, char(string(opts.YLabel)));
    end
    if strlength(string(opts.Title)) > 0
        title(ax, char(string(opts.Title)));
    end
end

function bNodes = interior_boundary_nodes(coords, bNodes)
    if isempty(bNodes) || size(coords, 2) < 2
        return;
    end

    x = coords(:, 1);
    y = coords(:, 2);
    xr = max(x) - min(x);
    yr = max(y) - min(y);
    tol = 1e-6 * max([xr, yr, 1]);

    onOuter = (abs(x - min(x)) <= tol) | (abs(x - max(x)) <= tol) | ...
              (abs(y - min(y)) <= tol) | (abs(y - max(y)) <= tol);
    bNodes = bNodes(~onOuter(bNodes));
end
