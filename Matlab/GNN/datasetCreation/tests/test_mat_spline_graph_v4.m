function test_mat_spline_graph_v4(varargin)
%TEST_MAT_SPLINE_GRAPH_V4 Regression checks for MAT/spline graph schema v4.
%
% Fast synthetic checks run by default:
%   test_mat_spline_graph_v4
%
% Optional real-data smoke test over three INP files:
%   test_mat_spline_graph_v4('RunRealSmoke', true)

    p = inputParser;
    p.addParameter('RunRealSmoke', false, @(x) islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    testsDir = fileparts(mfilename('fullpath'));
    datasetDir = fileparts(testsDir);
    gnnRoot = fileparts(datasetDir);
    repoRoot = fileparts(fileparts(gnnRoot));

    addpath(datasetDir);
    addpath(fullfile(datasetDir, 'gnn_graph'));
    addpath(genpath(fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal', 'src')));
    addpath(genpath(fullfile(gnnRoot, 'helpers')));

    run_synthetic_tests();
    if logical(opts.RunRealSmoke)
        run_real_smoke(repoRoot, gnnRoot, datasetDir);
    end

    fprintf('test_mat_spline_graph_v4 passed.\n');
end

function run_synthetic_tests()
    xVals = 1:64;
    yVals = 1:64;

    barMask = false(64, 64);
    barMask(30:34, 8:56) = true;
    check_case('bar', barMask, xVals, yVals, ...
        'MinEndpoints', 2, 'MinJunctions', 0);

    crossMask = false(64, 64);
    crossMask(29:35, 8:56) = true;
    crossMask(8:56, 29:35) = true;
    check_case('cross', crossMask, xVals, yVals, ...
        'MinEndpoints', 4, 'MinJunctions', 1);

    [xx, yy] = meshgrid(1:64, 1:64);
    rr = sqrt((xx - 32) .^ 2 + (yy - 32) .^ 2);
    loopMask = rr >= 12 & rr <= 19;
    check_case('loop', loopMask, xVals, yVals, ...
        'MinEndpoints', 0, 'MinJunctions', 0);

    islandMask = false(64, 64);
    islandMask(32, 32) = true;
    check_case('island', islandMask, xVals, yVals, ...
        'MinEndpoints', 0, 'MinJunctions', 0);
end

function check_case(name, mask, xVals, yVals, varargin)
    p = inputParser;
    p.addParameter('MinEndpoints', 0);
    p.addParameter('MinJunctions', 0);
    p.parse(varargin{:});
    opts = p.Results;

    [g, dbg] = mat_spline_graph_from_mask(mask, xVals, yVals);
    assert(g.num_nodes > 0, '%s: graph has no nodes.', name);
    assert(size(g.edge_index, 1) == 2, '%s: edge_index must be 2xE.', name);
    assert(size(g.node_features, 2) == 7, '%s: expected 7 v4 node features.', name);
    assert(all(isfinite(g.node_features(:, 1:4)), 'all'), ...
        '%s: node features contain non-finite geometry values.', name);
    assert(all(g.node_features(:, 3) >= 0), '%s: radius must be nonnegative.', name);
    assert(dbg.skeleton_component_count == dbg.occupancy_component_count, ...
        '%s: skeleton component count should preserve occupancy components.', name);
    assert(dbg.reconstruction_error <= dbg.error_tolerance + 1e-9, ...
        '%s: reconstruction error %.4g exceeds tolerance %.4g.', ...
        name, dbg.reconstruction_error, dbg.error_tolerance);

    nEndpoint = nnz(g.node_features(:, 5) > 0);
    nJunction = nnz(g.node_features(:, 6) > 0);
    assert(nEndpoint >= opts.MinEndpoints, ...
        '%s: expected at least %d endpoints, got %d.', name, opts.MinEndpoints, nEndpoint);
    assert(nJunction >= opts.MinJunctions, ...
        '%s: expected at least %d junctions, got %d.', name, opts.MinJunctions, nJunction);
end

function run_real_smoke(repoRoot, gnnRoot, datasetDir)
    samplesRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');
    inpFiles = dir(fullfile(samplesRoot, '**', 'sheet_shell.inp'));
    assert(numel(inpFiles) >= 3, 'Need at least 3 sheet_shell.inp files for real smoke test.');

    outRoot = tempname;
    samplesOut = fullfile(outRoot, 'samples');
    targetsOut = fullfile(outRoot, 'targets');
    mkdir(samplesOut);

    sampleIds = cell(3, 1);
    for i = 1:3
        inpPath = fullfile(inpFiles(i).folder, inpFiles(i).name);
        [~, runName] = fileparts(inpFiles(i).folder);
        sampleIds{i} = sprintf('%s_smoke%d', runName, i);

        inpData = read_abaqus_inp(inpPath);
        fullGraph = build_full_reference_graph_gnn(inpData, ...
            'ElsetName', 'SPINODAL_SHELL', ...
            'AutoDetectElset', true, ...
            'PatternPriority', {'spinodal', 'top'});

        [matGraph, matDebug] = extract_mat_spline_graph_gnn(fullGraph);
        [~, skeletonGraph] = extract_structural_graph_gnn(fullGraph, ...
            'DetailLevel', 1.00, ...
            'MinIslandNodes', 1);
        assert(matGraph.num_nodes < skeletonGraph.num_nodes, ...
            'v4 control graph should be smaller than detail-level-1 skeleton for %s.', inpPath);

        gnn_data = struct();
        gnn_data.edge_index = int32(matGraph.edge_index);
        gnn_data.num_nodes = matGraph.num_nodes;
        gnn_data.x = matGraph.node_features;
        gnn_data.feature_names = matGraph.feature_names;
        gnn_data.schema_version = 4;
        gnn_data.representation = 'mat_spline_control_graph';
        [gnn_data.tr_ratio, gnn_data.ang_deg] = parse_run_name(runName);

        dense_data = local_dense_raster(fullGraph, 128);
        matSplineGraph = matGraph; %#ok<NASGU>
        sampleDir = fullfile(samplesOut, sampleIds{i});
        mkdir(sampleDir);
        save(fullfile(sampleDir, 'sample.mat'), 'gnn_data', 'dense_data', 'matSplineGraph', 'matDebug', '-v7');
        loaded = load(fullfile(sampleDir, 'sample.mat'), 'gnn_data');
        assert(loaded.gnn_data.schema_version == 4, 'Reloaded sample schema mismatch.');
    end

    REPO_ROOT = repoRoot; %#ok<NASGU>
    SAMPLES_DIR = samplesOut; %#ok<NASGU>
    TARGETS_DIR = targetsOut; %#ok<NASGU>
    run(fullfile(datasetDir, 'step9_aggregate_mat_spline_graphs.m'));
    assert(isfile(fullfile(targetsOut, 'graphs_all.mat')), 'Smoke aggregate was not created.');

    [X_cell, ei_cell, N_vec, ~, Dense] = load_hybrid_graph_dataset(sampleIds, samplesOut);
    trainMask = true(numel(sampleIds), 1);
    [X_pad, nodeMask] = pad_and_normalize_hybrid_graphs(X_cell, N_vec, trainMask, []);
    A_hat = build_norm_adjacency(ei_cell, N_vec, max(N_vec));
    DenseInput = prepare_dense_raster_inputs(Dense);

    F = size(X_pad, 1);
    params = init_hybrid_params(F, 8, [4, 8, 16], 16, 2, 1, 3);
    globals = dlarray(zeros(3, 1, numel(sampleIds), 'single'));
    Zhat = hybrid_forward(params, X_pad, A_hat, 1, single(nodeMask), ...
        DenseInput, globals, false, 0, true, true);
    assert(isequal(size(Zhat), [2, 1, numel(sampleIds)]), 'Unexpected forward output size.');
end

function [tr_ratio, ang_deg] = parse_run_name(runName)
tk = regexp(runName, 'tr(\d+)', 'tokens', 'once');
tr_ratio = NaN;
if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end
tk = regexp(runName, 'ang(\d+)', 'tokens', 'once');
ang_deg = NaN;
if ~isempty(tk), ang_deg = str2double(tk{1}); end
end

function dense_data = local_dense_raster(fullGraph, gridSize)
xy = double(fullGraph.node_coords(:, 1:2));
boundary = logical(fullGraph.boundary_mask(:));
xLo = min(xy(:, 1)); xHi = max(xy(:, 1));
yLo = min(xy(:, 2)); yHi = max(xy(:, 2));
cols = 1 + round((xy(:, 1) - xLo) / max(xHi - xLo, eps) * (gridSize - 1));
rows = 1 + round((xy(:, 2) - yLo) / max(yHi - yLo, eps) * (gridSize - 1));
cols = min(gridSize, max(1, cols));
rows = min(gridSize, max(1, rows));
occupancy = false(gridSize, gridSize);
boundaryMask = false(gridSize, gridSize);
idx = sub2ind([gridSize, gridSize], rows, cols);
occupancy(idx) = true;
boundaryMask(idx(boundary)) = true;
dense_data = struct();
dense_data.raster = uint8(cat(3, occupancy, boundaryMask));
dense_data.channels = {'occupancy', 'boundary'};
dense_data.grid_size = gridSize;
dense_data.xy_bounds = [xLo, xHi, yLo, yHi];
dense_data.schema_version = 1;
end
