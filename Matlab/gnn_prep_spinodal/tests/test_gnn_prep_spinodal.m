function test_gnn_prep_spinodal()
%TEST_GNN_PREP_SPINODAL Minimal regression checks for gnn_prep_spinodal.

    testsDir = fileparts(mfilename('fullpath'));
    moduleRoot = fileparts(testsDir);
    matlabRoot = fileparts(moduleRoot);

    addpath(moduleRoot);
    addpath(genpath(fullfile(moduleRoot, 'src')));

    inpPath = fullfile(matlabRoot, 'results', 'sheets', 'lamellar', ...
        'sheetCone_tr50_ang000_lamellar_N128_1x1', 'FEA_shell', 'sheet_shell_shell_job.inp');
    csvPath = fullfile(matlabRoot, 'results', 'sheets', 'lamellar', ...
        'sheetCone_tr50_ang000_lamellar_N128_1x1', 'FEA_shell', 'midplane_results_shell.csv');

    assert(isfile(inpPath), 'Missing test INP: %s', inpPath);
    assert(isfile(csvPath), 'Missing test CSV: %s', csvPath);

    inpData = read_abaqus_inp(inpPath);
    fullGraph = build_full_reference_graph(inpData, 'ElsetName', 'SPINODAL_SHELL');
    fullGraph = attach_nodal_data(fullGraph, read_abaqus_nodal_csv(csvPath), 'Strict', false);

    assert(fullGraph.num_nodes > 0, 'Full graph has no nodes.');
    assert(size(fullGraph.edges_local, 1) > 0, 'Full graph has no edges.');
    assert(any(fullGraph.boundary_mask), 'Boundary detection produced no boundary nodes.');
    assert(height(fullGraph.node_data_table) == fullGraph.num_nodes, ...
        'Node data rows do not match full graph node count.');
    assert(isequal(fullGraph.node_data_table.Label, fullGraph.node_labels), ...
        'Label alignment mismatch after attach_nodal_data.');

    downA = downsample_graph_deterministic(fullGraph, 'DetailLevel', 0.20);
    downB = downsample_graph_deterministic(fullGraph, 'DetailLevel', 0.20);
    downGridA = downsample_graph_deterministic(fullGraph, ...
        'DetailLevel', 0.20, 'NodeSelectionMethod', 'grid');
    downGridB = downsample_graph_deterministic(fullGraph, ...
        'DetailLevel', 0.20, 'NodeSelectionMethod', 'grid');
    downHybridA = downsample_graph_deterministic(fullGraph, ...
        'DetailLevel', 0.20, 'NodeSelectionMethod', 'hybrid');
    downHybridB = downsample_graph_deterministic(fullGraph, ...
        'DetailLevel', 0.20, 'NodeSelectionMethod', 'hybrid');
    midExpA = downsample_midpoint_experimental(downHybridA);
    midExpB = downsample_midpoint_experimental(downHybridA);

    assert(downA.num_nodes <= fullGraph.num_nodes, 'Downsampled graph is larger than full graph.');
    assert(isequal(downA.node_labels, downB.node_labels), ...
        'Downsampling is not deterministic for identical inputs.');
    assert(isequal(downGridA.node_labels, downGridB.node_labels), ...
        'Grid downsampling is not deterministic for identical inputs.');
    assert(isequal(downHybridA.node_labels, downHybridB.node_labels), ...
        'Hybrid downsampling is not deterministic for identical inputs.');
    assert(isequal(midExpA.node_coords, midExpB.node_coords), ...
        'Experimental midpoint graph is not deterministic for identical inputs.');
    if midExpA.num_nodes > 0
        assert(isfield(midExpA, 'node_data_table') && any(strcmp(midExpA.node_data_table.Properties.VariableNames, 'VerticalGap')), ...
            'Experimental midpoint graph missing VerticalGap node attribute.');
        assert(isfield(midExpA, 'node_data_table') && any(strcmp(midExpA.node_data_table.Properties.VariableNames, 'Distance')), ...
            'Experimental midpoint graph missing Distance node attribute.');
        assert(all(midExpA.node_data_table.VerticalGap >= 0), ...
            'Experimental midpoint graph has negative VerticalGap.');
        if isfield(midExpA, 'midpoint_info')
            info = midExpA.midpoint_info;
            if isfield(info, 'n_midpoint_nodes_dense') && isfield(info, 'n_midpoint_nodes_sparse')
                assert(info.n_midpoint_nodes_sparse <= info.n_midpoint_nodes_dense, ...
                    'Experimental midpoint sparsification increased node count.');
            end
        end
    end
    assert(all(ismember(fullGraph.boundary_labels, downA.node_labels)), ...
        'Not all full-graph boundary nodes were preserved in downsample.');

    tmpOut = fullfile(moduleRoot, 'out', 'test_export');
    paths = export_graph_csv(downA, tmpOut, 'test_graph');
    assert(isfile(paths.nodes_csv_path), 'Missing exported nodes CSV.');
    assert(isfile(paths.edges_csv_path), 'Missing exported edges CSV.');
    assert(isfile(paths.mat_path), 'Missing exported MAT file.');

    fprintf('test_gnn_prep_spinodal passed.\n');
end
