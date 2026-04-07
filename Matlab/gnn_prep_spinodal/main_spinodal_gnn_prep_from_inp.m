function outputs = main_spinodal_gnn_prep_from_inp(inpPath, varargin)
%MAIN_SPINODAL_GNN_PREP_FROM_INP Build full/skeleton/structural graph from INP only.
%
% outputs = MAIN_SPINODAL_GNN_PREP_FROM_INP(inpPath, ...)
%
% This is the pre-deformation variant of main_spinodal_gnn_prep:
% - reads only Abaqus INP geometry
% - does not require nodal CSV

    if ~(ischar(inpPath) || isstring(inpPath))
        error('main_spinodal_gnn_prep_from_inp:BadInpPathType', ...
            'inpPath must be a char vector or string scalar.');
    end
    inpPath = char(string(inpPath));

    p = inputParser;
    p.addParameter('ElsetName', 'SPINODAL_SHELL', @(x) ischar(x) || isstring(x));
    p.addParameter('AutoDetectElset', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PatternPriority', {'spinodal', 'top'}, @(x) iscell(x) || isstring(x));
    p.addParameter('StructuralDetailLevel', 0.30, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('StructuralMinIslandNodes', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('OutDir', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Prefix', 'spinodal_graph_inp', @(x) ischar(x) || isstring(x));
    p.addParameter('RunName', '', @(x) ischar(x) || isstring(x));
    p.addParameter('MakePlots', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('SavePlots', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('PlotUnitsScale', 1000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    p.addParameter('PlotXLabel', '', @(x) ischar(x) || isstring(x));
    p.addParameter('PlotYLabel', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    opts = p.Results;

    moduleRoot = fileparts(mfilename('fullpath'));
    addpath(moduleRoot);
    addpath(genpath(fullfile(moduleRoot, 'src')));

    prefix = sanitize_token(char(string(opts.Prefix)));
    if isempty(prefix)
        prefix = 'spinodal_graph_inp';
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
    fullGraph = build_full_reference_graph(inpData, ...
        'ElsetName', opts.ElsetName, ...
        'AutoDetectElset', logical(opts.AutoDetectElset), ...
        'PatternPriority', opts.PatternPriority);

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
            'Title', 'Full Reference Graph (INP only)', ...
            'Style', 'paper2d', ...
            'BoundaryMode', 'interior', ...
            'UnitsScale', opts.PlotUnitsScale, ...
            'XLabel', xLabelText, ...
            'YLabel', yLabelText);

        figStruct = figure('Name', 'Structural Graph Overlay');
        axStruct = axes('Parent', figStruct);
        plot_structural_graph_overlay(fullGraph, skeletonGraph, structuralGraph, ...
            'AxesHandle', axStruct, ...
            'Title', sprintf('Structural Graph Overlay (detail=%.2f)', opts.StructuralDetailLevel), ...
            'UnitsScale', opts.PlotUnitsScale, ...
            'XLabel', xLabelText, ...
            'YLabel', yLabelText);

        if logical(opts.SavePlots)
            plotPaths.full_png = fullfile(dirs.full_plots, 'full_graph.png');
            plotPaths.structural_overlay_png = fullfile(dirs.struct_plots, 'structural_graph_overlay.png');
            saveas(figFull, plotPaths.full_png);
            saveas(figStruct, plotPaths.structural_overlay_png);
        end
    end

    write_run_metadata(fullfile(dirs.meta, 'run_info.txt'), inpPath, opts, ...
        fullGraph, skeletonGraph, structuralGraph, runDir);

    fprintf('[gnn_prep_spinodal][inp-only] Full graph: %d nodes, %d edges.\n', ...
        fullGraph.num_nodes, size(fullGraph.edges_local, 1));
    fprintf('[gnn_prep_spinodal][inp-only] Skeleton graph: %d nodes, %d edges.\n', ...
        skeletonGraph.num_nodes, size(skeletonGraph.edges_local, 1));
    fprintf('[gnn_prep_spinodal][inp-only] Structural graph: %d nodes, %d edges.\n', ...
        structuralGraph.num_nodes, size(structuralGraph.edges_local, 1));
    fprintf('[gnn_prep_spinodal][inp-only] Run directory: %s\n', runDir);

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

function write_run_metadata(metaPath, inpPath, opts, fullGraph, skeletonGraph, structuralGraph, runDir)
    fid = fopen(metaPath, 'w');
    if fid < 0
        warning('main_spinodal_gnn_prep_from_inp:MetaWriteFailed', ...
            'Could not write metadata file: %s', metaPath);
        return;
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'run_directory: %s\n', runDir);
    fprintf(fid, 'timestamp: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
    fprintf(fid, 'mode: inp_only_predeformation\n');
    fprintf(fid, 'inp_path: %s\n', inpPath);
    fprintf(fid, '\n[options]\n');
    fprintf(fid, 'elset_name: %s\n', char(string(opts.ElsetName)));
    fprintf(fid, 'structural_detail_level: %.6f\n', opts.StructuralDetailLevel);
    fprintf(fid, 'structural_min_island_nodes: %d\n', opts.StructuralMinIslandNodes);
    fprintf(fid, 'plot_units_scale: %.6g\n', opts.PlotUnitsScale);
    fprintf(fid, '\n[graph_sizes]\n');
    fprintf(fid, 'full_nodes: %d\n', fullGraph.num_nodes);
    fprintf(fid, 'full_edges: %d\n', size(fullGraph.edges_local, 1));
    fprintf(fid, 'skeleton_nodes: %d\n', skeletonGraph.num_nodes);
    fprintf(fid, 'skeleton_edges: %d\n', size(skeletonGraph.edges_local, 1));
    fprintf(fid, 'structural_nodes: %d\n', structuralGraph.num_nodes);
    fprintf(fid, 'structural_edges: %d\n', size(structuralGraph.edges_local, 1));
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
