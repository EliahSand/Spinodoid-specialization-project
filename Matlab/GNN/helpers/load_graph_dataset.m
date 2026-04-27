function [X_cell, ei_cell, N_vec] = load_graph_dataset(sample_ids, samplesDir)
%LOAD_GRAPH_DATASET Load structural graphs for a list of sample IDs.
%   Returns:
%     X_cell{g}  - 3xNg single  (features: x, y, radius)
%     ei_cell{g} - 2xEg double  (1-based edge indices)
%     N_vec(g)   - number of nodes in graph g
%
%   If data/dataset/targets/graphs_all.mat exists the aggregate is used
%   (loaded once per session via a persistent cache). Otherwise falls back
%   to loading individual sample.mat files.

persistent agg aggMap   % cached across calls in the same MATLAB session

G = numel(sample_ids);

%% Try aggregate path

aggFile = fullfile(fileparts(samplesDir), 'targets', 'graphs_all.mat');

if isfile(aggFile)

    % Load aggregate once per session
    if isempty(agg)
        fprintf('load_graph_dataset: loading aggregate from %s\n', aggFile);
        agg = load(aggFile, 'sample_ids', 'X_cell', 'ei_cell', 'N_vec');
        aggMap = containers.Map(agg.sample_ids, num2cell(1:numel(agg.sample_ids)));
    end

    X_cell  = cell(G, 1);
    ei_cell = cell(G, 1);
    N_vec   = zeros(G, 1);

    for g = 1:G
        sid = sample_ids{g};
        if ~isKey(aggMap, sid)
            error('load_graph_dataset:missing', ...
                'sample_id "%s" not found in graphs_all.mat. Re-run step9.', sid);
        end
        idx        = aggMap(sid);
        X_cell{g}  = agg.X_cell{idx}(1:3, :);   % drop boundary column (schema v2)
        ei_cell{g} = agg.ei_cell{idx};
        N_vec(g)   = agg.N_vec(idx);
    end

else
    % Fallback: load individual sample.mat files
    fprintf('load_graph_dataset: graphs_all.mat not found — loading %d files individually.\n', G);
    fprintf('  (run step9_aggregate_graphs.m to speed this up)\n');

    X_cell  = cell(G, 1);
    ei_cell = cell(G, 1);
    N_vec   = zeros(G, 1, 'int32');

    for g = 1:G
        sid  = sample_ids{g};
        path = fullfile(samplesDir, sid, 'sample.mat');
        if ~isfile(path)
            error('load_graph_dataset:missing', 'sample.mat not found: %s', path);
        end
        d  = load(path, 'gnn_data');
        gd = d.gnn_data;

        X_cell{g}  = single(gd.x(:, 1:3).');   % drop boundary column (schema v2)
        ei_cell{g} = double(gd.edge_index);
        N_vec(g)   = gd.num_nodes;
    end

    N_vec = double(N_vec);
end
end
