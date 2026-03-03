function outputs = main_spinodal_gnn_prep(inpPath, csvPath, varargin)
%MAIN_SPINODAL_GNN_PREP End-to-end Abaqus-to-graph preprocessing pipeline.
%
% outputs = MAIN_SPINODAL_GNN_PREP(inpPath, csvPath, ...)
%
% Name-Value options:
%   'ElsetName'            : explicit elset for region selection (default '')
%   'AutoDetectElset'      : auto-detect spinodal/top elset (default true)
%   'PatternPriority'      : auto-detect keyword priority
%   'DetailLevel'          : [0,1] downsampling detail level (default 0.25)
%   'TargetNumNodes'       : explicit downsample target nodes (optional)
%   'PreserveConnectivity' : shortest-path augmentation (default true)
%   'DownsampleMethod'     : 'farthest'|'grid'|'hybrid' (default 'farthest')
%   'HybridAnchorFraction' : [0,1], anchor fraction in hybrid mode (default 0.35)
%   'BoundaryKeepRatio'    : [0,1], keep ratio for boundary nodes (default 1)
%   'BoundaryStraightTolDeg' : corner tolerance for boundary simplification (default 8)
%   'OutDir'               : output directory (default ./out)
%   'Prefix'               : filename prefix (default 'spinodal_graph')
%   'RunName'              : explicit run folder name (optional)
%   'MakePlots'            : generate full/downsample plots (default true)
%   'SavePlots'            : write PNGs to OutDir (default true)
%   'StrictCSVMapping'     : error on missing CSV rows (default false)
%   'DownsampleOutputMode' : 'mesh' or 'line' (default 'mesh')
%   'LineNodeBudget'       : node budget in line mode (optional)
%   'PlotStyle'            : 'default' or 'paper2d' (default 'paper2d')
%   'PlotBoundaryMode'     : 'all'|'interior'|'none' (default 'interior')
%   'PlotUnitsScale'       : axis scale factor (default 1000; m->mm)
%   'PlotXLabel'           : explicit x-axis label override (optional)
%   'PlotYLabel'           : explicit y-axis label override (optional)
%   'PlotNodeMarkerSize'     : regular node marker size (optional)
%   'PlotBoundaryMarkerSize' : red boundary marker size (optional)

    if ~(ischar(inpPath) || isstring(inpPath))
        error('main_spinodal_gnn_prep:BadInpPathType', ...
            'inpPath must be a char vector or string scalar.');
    end
    if ~(ischar(csvPath) || isstring(csvPath))
        error('main_spinodal_gnn_prep:BadCsvPathType', ...
            'csvPath must be a char vector or string scalar.');
    end
    inpPath = char(string(inpPath));
    csvPath = char(string(csvPath));

    p = inputParser;
    p.addParameter('ElsetName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoDetectElset', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PatternPriority', {'spinodal', 'top'}, @(x) iscell(x) || isstring(x));
    p.addParameter('DetailLevel', 0.25, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('TargetNumNodes', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    p.addParameter('PreserveConnectivity', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('DownsampleMethod', 'farthest', @(x) ischar(x) || isstring(x));
    p.addParameter('HybridAnchorFraction', 0.35, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('BoundaryKeepRatio', 1, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.addParameter('BoundaryStraightTolDeg', 8, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 180);
    p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Prefix', 'spinodal_graph', @(x) ischar(x) || isstring(x));
    p.addParameter('RunName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('MakePlots', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('SavePlots', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('StrictCSVMapping', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('DownsampleOutputMode', 'mesh', @(x) ischar(x) || isstring(x));
    p.addParameter('LineNodeBudget', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    p.addParameter('PlotStyle', 'paper2d', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotBoundaryMode', 'interior', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotUnitsScale', 1000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('PlotXLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotYLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotNodeMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.addParameter('PlotBoundaryMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.parse(varargin{:});
    opts = p.Results;
    downOutputMode = lower(char(string(opts.DownsampleOutputMode)));
    downOutputMode = validatestring(downOutputMode, {'mesh', 'line'});
    downsampleMethod = lower(char(string(opts.DownsampleMethod)));
    downsampleMethod = validatestring(downsampleMethod, {'farthest', 'grid', 'hybrid'});
    opts.DownsampleMethod = downsampleMethod;

    moduleRoot = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(moduleRoot, 'src')));

    prefix = sanitize_token(char(string(opts.Prefix)));
    if isempty(prefix)
        prefix = 'spinodal_graph';
    end

    if strlength(string(opts.OutDir)) == 0
        outRoot = fullfile(moduleRoot, 'out');
    else
        outRoot = char(string(opts.OutDir));
    end
    ensure_dir(outRoot);

    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    if strlength(string(opts.RunName)) > 0
        runFolder = sanitize_token(char(string(opts.RunName)));
        if isempty(runFolder)
            runFolder = sprintf('%s_%s', prefix, stamp);
        end
    else
        runFolder = sprintf('%s_%s', prefix, stamp);
    end
    runDir = fullfile(outRoot, runFolder);

    dirs = struct();
    dirs.run = runDir;
    dirs.meta = fullfile(runDir, 'metadata');
    dirs.full_data = fullfile(runDir, 'full_graph', 'data');
    dirs.full_plots = fullfile(runDir, 'full_graph', 'plots');
    dirs.down_data = fullfile(runDir, 'downsampled_graph', 'data');
    dirs.down_plots = fullfile(runDir, 'downsampled_graph', 'plots');
    dirs.line_data = fullfile(runDir, 'line_graph', 'data');
    dirs.line_plots = fullfile(runDir, 'line_graph', 'plots');
    ensure_dir(dirs.meta);
    ensure_dir(dirs.full_data);
    ensure_dir(dirs.full_plots);
    ensure_dir(dirs.down_data);
    ensure_dir(dirs.down_plots);
    if strcmp(downOutputMode, 'line')
        ensure_dir(dirs.line_data);
        ensure_dir(dirs.line_plots);
    end

    inpData = read_abaqus_inp(inpPath);
    nodalTable = read_abaqus_nodal_csv(csvPath);

    fullGraph = build_full_reference_graph(inpData, ...
        'ElsetName', opts.ElsetName, ...
        'AutoDetectElset', logical(opts.AutoDetectElset), ...
        'PatternPriority', opts.PatternPriority);

    fullGraph = attach_nodal_data(fullGraph, nodalTable, ...
        'Strict', logical(opts.StrictCSVMapping));

    downMeshGraph = downsample_graph_deterministic(fullGraph, ...
        'DetailLevel', opts.DetailLevel, ...
        'TargetNumNodes', opts.TargetNumNodes, ...
        'PreserveConnectivity', logical(opts.PreserveConnectivity), ...
        'NodeSelectionMethod', opts.DownsampleMethod, ...
        'HybridAnchorFraction', opts.HybridAnchorFraction, ...
        'BoundaryKeepRatio', opts.BoundaryKeepRatio, ...
        'BoundaryStraightTolDeg', opts.BoundaryStraightTolDeg, ...
        'OutputMode', 'mesh');

    lineGraph = struct();
    downGraph = downMeshGraph;
    if strcmp(downOutputMode, 'line')
        lineGraph = downsample_graph_deterministic(fullGraph, ...
            'DetailLevel', opts.DetailLevel, ...
            'TargetNumNodes', opts.TargetNumNodes, ...
            'PreserveConnectivity', logical(opts.PreserveConnectivity), ...
            'NodeSelectionMethod', opts.DownsampleMethod, ...
            'HybridAnchorFraction', opts.HybridAnchorFraction, ...
            'BoundaryKeepRatio', opts.BoundaryKeepRatio, ...
            'BoundaryStraightTolDeg', opts.BoundaryStraightTolDeg, ...
            'OutputMode', 'line', ...
            'LineNodeBudget', opts.LineNodeBudget);
        downGraph = lineGraph;
    end

    fullPaths = export_graph_csv(fullGraph, dirs.full_data, 'full_graph');
    downPaths = export_graph_csv(downMeshGraph, dirs.down_data, 'downsampled_graph');
    linePaths = struct();
    if strcmp(downOutputMode, 'line')
        linePaths = export_graph_csv(lineGraph, dirs.line_data, 'line_graph');
    end

    xLabelText = char(string(opts.PlotXLabel));
    yLabelText = char(string(opts.PlotYLabel));
    if isempty(strtrim(xLabelText))
        xLabelText = x_axis_label(opts.PlotUnitsScale);
    end
    if isempty(strtrim(yLabelText))
        yLabelText = y_axis_label(opts.PlotUnitsScale);
    end

    plotPaths = struct();
    if logical(opts.MakePlots)
        figFull = figure('Name', 'Full Reference Graph');
        axFull = axes('Parent', figFull);
        plot_graph_mesh(fullGraph, ...
            'AxesHandle', axFull, ...
            'Title', 'Full Reference Graph', ...
            'Style', opts.PlotStyle, ...
            'BoundaryMode', opts.PlotBoundaryMode, ...
            'UnitsScale', opts.PlotUnitsScale, ...
            'XLabel', xLabelText, ...
            'YLabel', yLabelText, ...
            'NodeMarkerSize', opts.PlotNodeMarkerSize, ...
            'BoundaryMarkerSize', opts.PlotBoundaryMarkerSize);

        figDown = figure('Name', 'Downsampled Graph');
        axDown = axes('Parent', figDown);
        plot_graph_mesh(downMeshGraph, ...
            'AxesHandle', axDown, ...
            'Title', sprintf('Downsampled Graph (detail=%.2f)', opts.DetailLevel), ...
            'Style', opts.PlotStyle, ...
            'BoundaryMode', opts.PlotBoundaryMode, ...
            'UnitsScale', opts.PlotUnitsScale, ...
            'XLabel', xLabelText, ...
            'YLabel', yLabelText, ...
            'NodeMarkerSize', opts.PlotNodeMarkerSize, ...
            'BoundaryMarkerSize', opts.PlotBoundaryMarkerSize);

        figLine = [];
        if strcmp(downOutputMode, 'line')
            figLine = figure('Name', 'Line Graph');
            axLine = axes('Parent', figLine);
            plot_graph_mesh(lineGraph, ...
                'AxesHandle', axLine, ...
                'Title', 'Line Graph', ...
                'Style', opts.PlotStyle, ...
                'BoundaryMode', 'none', ...
                'UnitsScale', opts.PlotUnitsScale, ...
                'XLabel', xLabelText, ...
                'YLabel', yLabelText, ...
                'NodeMarkerSize', opts.PlotNodeMarkerSize, ...
                'BoundaryMarkerSize', opts.PlotBoundaryMarkerSize);
        end

        if logical(opts.SavePlots)
            plotPaths.full_png = fullfile(dirs.full_plots, 'full_graph.png');
            plotPaths.down_png = fullfile(dirs.down_plots, 'downsampled_graph.png');
            saveas(figFull, plotPaths.full_png);
            saveas(figDown, plotPaths.down_png);
            if ~isempty(figLine)
                plotPaths.line_png = fullfile(dirs.line_plots, 'line_graph.png');
                saveas(figLine, plotPaths.line_png);
            end
        end
    end

    write_run_metadata(fullfile(dirs.meta, 'run_info.txt'), inpPath, csvPath, opts, ...
        fullGraph, downMeshGraph, downGraph, runDir);

    fprintf('[gnn_prep_spinodal] Full graph: %d nodes, %d edges, %d boundary nodes.\n', ...
        fullGraph.num_nodes, size(fullGraph.edges_local, 1), nnz(fullGraph.boundary_mask));
    fprintf('[gnn_prep_spinodal] Downsampled mesh graph: %d nodes, %d edges.\n', ...
        downMeshGraph.num_nodes, size(downMeshGraph.edges_local, 1));
    if strcmp(downOutputMode, 'line')
        fprintf('[gnn_prep_spinodal] Line graph: %d nodes, %d edges.\n', ...
            downGraph.num_nodes, size(downGraph.edges_local, 1));
    end
    fprintf('[gnn_prep_spinodal] Run directory: %s\n', runDir);

    outputs = struct();
    outputs.run_dir = runDir;
    outputs.out_dirs = dirs;
    outputs.inp_data = inpData;
    outputs.full_graph = fullGraph;
    outputs.downsampled_mesh_graph = downMeshGraph;
    if strcmp(downOutputMode, 'line')
        outputs.line_graph = downGraph;
    end
    outputs.down_graph = downGraph;
    outputs.exports = struct('full', fullPaths, 'down', downPaths, 'line', linePaths, 'plots', plotPaths);
end

function ensure_dir(pathStr)
    if ~isfolder(pathStr)
        mkdir(pathStr);
    end
end

function out = sanitize_token(in)
    out = regexprep(in, '[^A-Za-z0-9_-]+', '_');
    out = regexprep(out, '_+', '_');
    out = strtrim(out);
end

function write_run_metadata(metaPath, inpPath, csvPath, opts, fullGraph, downMeshGraph, downGraph, runDir)
    fid = fopen(metaPath, 'w');
    if fid < 0
        warning('main_spinodal_gnn_prep:MetaWriteFailed', 'Could not write metadata file: %s', metaPath);
        return;
    end
    c = onCleanup(@() fclose(fid));

    fprintf(fid, 'run_directory: %s\n', runDir);
    fprintf(fid, 'timestamp: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'inp_path: %s\n', inpPath);
    fprintf(fid, 'csv_path: %s\n', csvPath);
    fprintf(fid, '\n[options]\n');
    fprintf(fid, 'elset_name: %s\n', char(string(opts.ElsetName)));
    fprintf(fid, 'detail_level: %.6f\n', opts.DetailLevel);
    fprintf(fid, 'target_num_nodes: %s\n', mat2str(opts.TargetNumNodes));
    fprintf(fid, 'preserve_connectivity: %d\n', logical(opts.PreserveConnectivity));
    fprintf(fid, 'downsample_method: %s\n', char(string(opts.DownsampleMethod)));
    fprintf(fid, 'hybrid_anchor_fraction: %.6f\n', opts.HybridAnchorFraction);
    fprintf(fid, 'boundary_keep_ratio: %.6f\n', opts.BoundaryKeepRatio);
    fprintf(fid, 'boundary_straight_tol_deg: %.6f\n', opts.BoundaryStraightTolDeg);
    fprintf(fid, 'downsample_output_mode: %s\n', char(string(opts.DownsampleOutputMode)));
    fprintf(fid, 'line_node_budget: %s\n', mat2str(opts.LineNodeBudget));
    fprintf(fid, 'plot_style: %s\n', char(string(opts.PlotStyle)));
    fprintf(fid, 'plot_boundary_mode: %s\n', char(string(opts.PlotBoundaryMode)));
    fprintf(fid, 'plot_units_scale: %.6g\n', opts.PlotUnitsScale);
    fprintf(fid, 'plot_x_label: %s\n', char(string(opts.PlotXLabel)));
    fprintf(fid, 'plot_y_label: %s\n', char(string(opts.PlotYLabel)));
    fprintf(fid, 'plot_node_marker_size: %s\n', mat2str(opts.PlotNodeMarkerSize));
    fprintf(fid, 'plot_boundary_marker_size: %s\n', mat2str(opts.PlotBoundaryMarkerSize));

    fprintf(fid, '\n[graph_sizes]\n');
    fprintf(fid, 'full_nodes: %d\n', fullGraph.num_nodes);
    fprintf(fid, 'full_edges: %d\n', size(fullGraph.edges_local, 1));
    fprintf(fid, 'downsampled_mesh_nodes: %d\n', downMeshGraph.num_nodes);
    fprintf(fid, 'downsampled_mesh_edges: %d\n', size(downMeshGraph.edges_local, 1));
    fprintf(fid, 'down_output_nodes: %d\n', downGraph.num_nodes);
    fprintf(fid, 'down_output_edges: %d\n', size(downGraph.edges_local, 1));
end

function txt = x_axis_label(unitsScale)
    txt = axis_label_with_units('x', unitsScale);
end

function txt = y_axis_label(unitsScale)
    txt = axis_label_with_units('y', unitsScale);
end

function txt = axis_label_with_units(axisName, unitsScale)
    if abs(unitsScale - 1000) < 1e-12
        txt = sprintf('%s [mm]', axisName);
    elseif abs(unitsScale - 1) < 1e-12
        txt = sprintf('%s [m]', axisName);
    else
        txt = axisName;
    end
end
