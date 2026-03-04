function outputs = main_spinodal_gnn_prep(inpPath, csvPath, varargin)
%MAIN_SPINODAL_GNN_PREP End-to-end Abaqus-to-graph preprocessing pipeline.
%
% outputs = MAIN_SPINODAL_GNN_PREP(inpPath, csvPath, ...)
%
% This streamlined pipeline keeps:
%   1. the full reference graph
%   2. the extracted skeleton graph
%   3. the reduced structural graph
%
% Name-Value options:
%   'ElsetName'               : explicit elset for region selection (default '')
%   'AutoDetectElset'         : auto-detect spinodal/top elset (default true)
%   'PatternPriority'         : auto-detect keyword priority
%   'StructuralDetailLevel'   : [0,1], reduced structural graph detail (default 0.20)
%   'StructuralMinIslandNodes': keep at least this many skeleton nodes per component (default 1)
%   'OutDir'                  : output directory (default ./out)
%   'Prefix'                  : filename prefix (default 'spinodal_graph')
%   'RunName'                 : explicit run folder name (optional)
%   'MakePlots'               : generate full/structural plots (default true)
%   'SavePlots'               : write PNGs to OutDir (default true)
%   'StrictCSVMapping'        : error on missing CSV rows (default false)
%   'PlotStyle'               : 'default' or 'paper2d' (default 'paper2d')
%   'PlotBoundaryMode'        : 'all'|'interior'|'none' (default 'interior')
%   'MakeThicknessPlot'       : make separate dense+structural thickness plot (default false)
%   'PlotUnitsScale'          : axis scale factor (default 1000; m->mm)
%   'PlotXLabel'              : explicit x-axis label override (optional)
%   'PlotYLabel'              : explicit y-axis label override (optional)
%   'PlotNodeMarkerSize'      : regular node marker size (optional)
%   'PlotBoundaryMarkerSize'  : red boundary marker size (optional)

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
    p.addParameter('StructuralDetailLevel', 0.20, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('StructuralMinIslandNodes', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Prefix', 'spinodal_graph', @(x) ischar(x) || isstring(x));
    p.addParameter('RunName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('MakePlots', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('SavePlots', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('StrictCSVMapping', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PlotStyle', 'paper2d', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotBoundaryMode', 'interior', @(x) ischar(x) || isstring(x));
    p.addParameter('MakeThicknessPlot', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PlotUnitsScale', 1000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('PlotXLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotYLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotNodeMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.addParameter('PlotBoundaryMarkerSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    p.parse(varargin{:});
    opts = p.Results;

    moduleRoot = fileparts(mfilename('fullpath'));
    addpath(moduleRoot);
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
    dirs.struct_data = fullfile(runDir, 'structural_graph', 'data');
    dirs.struct_plots = fullfile(runDir, 'structural_graph', 'plots');
    ensure_dir(dirs.meta);
    ensure_dir(dirs.full_data);
    ensure_dir(dirs.full_plots);
    ensure_dir(dirs.struct_data);
    ensure_dir(dirs.struct_plots);

    inpData = read_abaqus_inp(inpPath);
    nodalTable = read_abaqus_nodal_csv(csvPath);

    fullGraph = build_full_reference_graph(inpData, ...
        'ElsetName', opts.ElsetName, ...
        'AutoDetectElset', logical(opts.AutoDetectElset), ...
        'PatternPriority', opts.PatternPriority);

    fullGraph = attach_nodal_data(fullGraph, nodalTable, ...
        'Strict', logical(opts.StrictCSVMapping));

    [structuralGraph, skeletonGraph, structuralDebug, fullGraphAnnotated] = extract_structural_graph(fullGraph, ...
        'DetailLevel', opts.StructuralDetailLevel, ...
        'MinIslandNodes', opts.StructuralMinIslandNodes);
    fullGraph = fullGraphAnnotated;

    fullPaths = export_graph_csv(fullGraph, dirs.full_data, 'full_graph');
    structuralPaths = struct();
    structuralPaths.skeleton = export_graph_csv(skeletonGraph, dirs.struct_data, 'skeleton_graph');
    structuralPaths.structural = export_graph_csv(structuralGraph, dirs.struct_data, 'structural_graph');

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

        figStruct = figure('Name', 'Structural Graph Overlay');
        axStruct = axes('Parent', figStruct);
        plot_structural_graph_overlay(fullGraph, skeletonGraph, structuralGraph, ...
            'AxesHandle', axStruct, ...
            'Title', sprintf('Structural Graph Overlay (detail=%.2f)', opts.StructuralDetailLevel), ...
            'UnitsScale', opts.PlotUnitsScale, ...
            'XLabel', xLabelText, ...
            'YLabel', yLabelText);

        figThickness = [];
        if logical(opts.MakeThicknessPlot)
            figThickness = figure('Name', 'Structural Thickness Overlay');
            axThickness = axes('Parent', figThickness);
            plot_structural_thickness_overlay(fullGraph, structuralGraph, ...
                'AxesHandle', axThickness, ...
                'Title', 'Structural Thickness Overlay', ...
                'UnitsScale', opts.PlotUnitsScale, ...
                'XLabel', xLabelText, ...
                'YLabel', yLabelText);
        end

        if logical(opts.SavePlots)
            plotPaths.full_png = fullfile(dirs.full_plots, 'full_graph.png');
            plotPaths.structural_overlay_png = fullfile(dirs.struct_plots, 'structural_graph_overlay.png');
            saveas(figFull, plotPaths.full_png);
            saveas(figStruct, plotPaths.structural_overlay_png);
            if ~isempty(figThickness)
                plotPaths.structural_thickness_overlay_png = fullfile(dirs.struct_plots, 'structural_thickness_overlay.png');
                saveas(figThickness, plotPaths.structural_thickness_overlay_png);
            end
        end
    end

    write_run_metadata(fullfile(dirs.meta, 'run_info.txt'), inpPath, csvPath, opts, ...
        fullGraph, skeletonGraph, structuralGraph, runDir);

    fprintf('[gnn_prep_spinodal] Full graph: %d nodes, %d edges, %d boundary nodes.\n', ...
        fullGraph.num_nodes, size(fullGraph.edges_local, 1), nnz(fullGraph.boundary_mask));
    fprintf('[gnn_prep_spinodal] Skeleton graph: %d nodes, %d edges.\n', ...
        skeletonGraph.num_nodes, size(skeletonGraph.edges_local, 1));
    fprintf('[gnn_prep_spinodal] Structural graph: %d nodes, %d edges.\n', ...
        structuralGraph.num_nodes, size(structuralGraph.edges_local, 1));
    fprintf('[gnn_prep_spinodal] Run directory: %s\n', runDir);

    outputs = struct();
    outputs.run_dir = runDir;
    outputs.out_dirs = dirs;
    outputs.inp_data = inpData;
    outputs.full_graph = fullGraph;
    outputs.skeleton_graph = skeletonGraph;
    outputs.structural_graph = structuralGraph;
    outputs.structural_debug = structuralDebug;
    outputs.exports = struct( ...
        'full', fullPaths, ...
        'structural', structuralPaths, ...
        'plots', plotPaths);
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

function write_run_metadata(metaPath, inpPath, csvPath, opts, fullGraph, skeletonGraph, structuralGraph, runDir)
    fid = fopen(metaPath, 'w');
    if fid < 0
        warning('main_spinodal_gnn_prep:MetaWriteFailed', 'Could not write metadata file: %s', metaPath);
        return;
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'run_directory: %s\n', runDir);
    fprintf(fid, 'timestamp: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'inp_path: %s\n', inpPath);
    fprintf(fid, 'csv_path: %s\n', csvPath);
    fprintf(fid, '\n[options]\n');
    fprintf(fid, 'elset_name: %s\n', char(string(opts.ElsetName)));
    fprintf(fid, 'structural_detail_level: %.6f\n', opts.StructuralDetailLevel);
    fprintf(fid, 'structural_min_island_nodes: %d\n', opts.StructuralMinIslandNodes);
    fprintf(fid, 'plot_style: %s\n', char(string(opts.PlotStyle)));
    fprintf(fid, 'plot_boundary_mode: %s\n', char(string(opts.PlotBoundaryMode)));
    fprintf(fid, 'make_thickness_plot: %d\n', logical(opts.MakeThicknessPlot));
    fprintf(fid, 'plot_units_scale: %.6g\n', opts.PlotUnitsScale);
    fprintf(fid, 'plot_x_label: %s\n', char(string(opts.PlotXLabel)));
    fprintf(fid, 'plot_y_label: %s\n', char(string(opts.PlotYLabel)));
    fprintf(fid, 'plot_node_marker_size: %s\n', mat2str(opts.PlotNodeMarkerSize));
    fprintf(fid, 'plot_boundary_marker_size: %s\n', mat2str(opts.PlotBoundaryMarkerSize));

    fprintf(fid, '\n[graph_sizes]\n');
    fprintf(fid, 'full_nodes: %d\n', fullGraph.num_nodes);
    fprintf(fid, 'full_edges: %d\n', size(fullGraph.edges_local, 1));
    fprintf(fid, 'skeleton_nodes: %d\n', skeletonGraph.num_nodes);
    fprintf(fid, 'skeleton_edges: %d\n', size(skeletonGraph.edges_local, 1));
    fprintf(fid, 'structural_nodes: %d\n', structuralGraph.num_nodes);
    fprintf(fid, 'structural_edges: %d\n', size(structuralGraph.edges_local, 1));
    if isfield(structuralGraph, 'cleanup_info') && isstruct(structuralGraph.cleanup_info)
        fprintf(fid, '\n[cleanup]\n');
        fprintf(fid, 'grid_spacing: %.6g\n', structuralGraph.cleanup_info.grid_spacing);
        fprintf(fid, 'min_edge_length: %.6g\n', structuralGraph.cleanup_info.min_edge_length);
        fprintf(fid, 'merge_radius: %.6g\n', structuralGraph.cleanup_info.merge_radius);
        fprintf(fid, 'iterations_used: %d\n', structuralGraph.cleanup_info.iterations_used);
        fprintf(fid, 'short_edge_contractions: %d\n', structuralGraph.cleanup_info.short_edge_contractions);
        fprintf(fid, 'cluster_merges: %d\n', structuralGraph.cleanup_info.cluster_merges);
    end
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
