function test_gnn_prep_spinodal()
%TEST_GNN_PREP_SPINODAL Minimal regression checks for the structural pipeline.

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
    [structGraphA, skeletonGraphA, debugA, fullGraphA] = extract_structural_graph(fullGraph, 'DetailLevel', 0.20);
    [structGraphB, skeletonGraphB] = extract_structural_graph(fullGraph, 'DetailLevel', 0.20);
    fullGraph = fullGraphA;

    assert(fullGraph.num_nodes > 0, 'Full graph has no nodes.');
    assert(size(fullGraph.edges_local, 1) > 0, 'Full graph has no edges.');
    assert(any(fullGraph.boundary_mask), 'Boundary detection produced no boundary nodes.');
    assert(height(fullGraph.node_data_table) == fullGraph.num_nodes, ...
        'Node data rows do not match full graph node count.');
    assert(isequal(fullGraph.node_data_table.Label, fullGraph.node_labels), ...
        'Label alignment mismatch after attach_nodal_data.');
    assert(isfield(fullGraph, 'node_thickness') && numel(fullGraph.node_thickness) == fullGraph.num_nodes, ...
        'Full graph is missing node thickness.');
    assert(isfield(fullGraph, 'node_radius') && numel(fullGraph.node_radius) == fullGraph.num_nodes, ...
        'Full graph is missing node radius.');

    assert(skeletonGraphA.num_nodes > 0, 'Skeleton graph has no nodes.');
    assert(structGraphA.num_nodes > 0, 'Structural graph has no nodes.');
    assert(structGraphA.num_nodes <= skeletonGraphA.num_nodes, ...
        'Structural graph should not have more nodes than skeleton graph.');
    assert(size(structGraphA.edge_index, 1) == 2, ...
        'Structural graph edge_index must be 2xE.');
    assert(numel(structGraphA.reduced_node_to_full_indices) == structGraphA.num_nodes, ...
        'Structural graph node-to-full mapping size mismatch.');
    assert(numel(structGraphA.edge_polyline) == size(structGraphA.edges_local, 1), ...
        'Structural graph edge polyline count mismatch.');
    assert(isfield(structGraphA, 'cleanup_info') && isstruct(structGraphA.cleanup_info), ...
        'Structural graph is missing cleanup metadata.');
    assert(isfield(skeletonGraphA, 'node_thickness') && numel(skeletonGraphA.node_thickness) == skeletonGraphA.num_nodes, ...
        'Skeleton graph is missing node thickness.');
    assert(isfield(structGraphA, 'node_thickness') && numel(structGraphA.node_thickness) == structGraphA.num_nodes, ...
        'Structural graph is missing node thickness.');
    assert(isfield(structGraphA, 'edge_thickness_min') && numel(structGraphA.edge_thickness_min) == size(structGraphA.edges_local, 1), ...
        'Structural graph is missing edge thickness summaries.');
    assert(isfield(debugA, 'structural_before_cleanup'), ...
        'Structural debug output is missing pre-cleanup graph.');
    assert(structGraphA.num_nodes <= debugA.structural_before_cleanup.num_nodes, ...
        'Cleanup should not increase structural node count.');
    if ~isempty(structGraphA.edge_length)
        assert(all(structGraphA.edge_length >= 0), ...
            'Structural graph edge lengths must be non-negative.');
    end

    assert(isequal(skeletonGraphA.node_labels, skeletonGraphB.node_labels), ...
        'Skeleton extraction is not deterministic for identical inputs.');
    assert(isequal(structGraphA.node_labels, structGraphB.node_labels), ...
        'Structural extraction is not deterministic for identical inputs.');
    assert(isequal(structGraphA.edges_local, structGraphB.edges_local), ...
        'Structural graph edges are not deterministic for identical inputs.');

    tmpOut = fullfile(moduleRoot, 'out', 'test_export');
    fullPaths = export_graph_csv(fullGraph, tmpOut, 'test_full_graph');
    skeletonPaths = export_graph_csv(skeletonGraphA, tmpOut, 'test_skeleton_graph');
    structPaths = export_graph_csv(structGraphA, tmpOut, 'test_structural_graph');
    assert(isfile(fullPaths.nodes_csv_path), 'Missing exported full nodes CSV.');
    assert(isfile(fullPaths.edges_csv_path), 'Missing exported full edges CSV.');
    assert(isfile(fullPaths.mat_path), 'Missing exported full MAT file.');
    assert(isfile(skeletonPaths.nodes_csv_path), 'Missing exported skeleton nodes CSV.');
    assert(isfile(skeletonPaths.edges_csv_path), 'Missing exported skeleton edges CSV.');
    assert(isfile(skeletonPaths.mat_path), 'Missing exported skeleton MAT file.');
    assert(isfile(structPaths.nodes_csv_path), 'Missing exported structural nodes CSV.');
    assert(isfile(structPaths.edges_csv_path), 'Missing exported structural edges CSV.');
    assert(isfile(structPaths.mat_path), 'Missing exported structural MAT file.');

    outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
        'ElsetName', 'SPINODAL_SHELL', ...
        'StructuralDetailLevel', 0.20, ...
        'MakePlots', false, ...
        'SavePlots', false, ...
        'OutDir', fullfile(moduleRoot, 'out'), ...
        'Prefix', 'test_pipeline');
    assert(isfield(outputs, 'full_graph') && isfield(outputs, 'skeleton_graph') && isfield(outputs, 'structural_graph'), ...
        'Main pipeline outputs are missing required graph products.');

    fprintf('test_gnn_prep_spinodal passed.\n');
end
