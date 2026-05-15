function [X_cell, ei_cell, N_vec, G_mat, Dense, EA_cell, datasetVersion] = load_hybrid_graph_dataset(sample_ids, samplesDir)
%LOAD_HYBRID_GRAPH_DATASET Load graph + dense raster inputs for hybrid models.
%   X_cell{g}        - F x Ng single. V1: [x; y; radius; boundary]   (F=4)
%                                     V2: [x; y; radius; bc_oneHot4] (F=7)
%   ei_cell{g}       - 2 x Eg double edge index
%   N_vec(g)         - node count
%   G_mat            - G x 2 [tr_ratio, ang_deg]
%   Dense            - H x W x 2 x G uint8 [occupancy, boundary]
%   EA_cell{g}       - 5 x Eg single edge_attr (V2). Empty for V1 graphs.
%   datasetVersion   - 'V1' or 'V2' character vector

persistent agg aggMap hasDense hasGlobals hasEdgeAttr versionLabel

G = numel(sample_ids);
aggFile = fullfile(fileparts(samplesDir), 'targets', 'graphs_all.mat');

if isfile(aggFile)
    if isempty(agg)
        fprintf('load_hybrid_graph_dataset: loading aggregate from %s\n', aggFile);
        vars = who('-file', aggFile);
        hasDense    = ismember('Dense_cell', vars);
        hasGlobals  = ismember('TR_vec', vars) && ismember('ANG_vec', vars);
        hasEdgeAttr = ismember('EA_cell', vars);
        hasVerLabel = ismember('schema_version_label', vars);
        loadVars = {'sample_ids', 'X_cell', 'ei_cell', 'N_vec'};
        if hasDense,    loadVars{end+1} = 'Dense_cell'; end %#ok<AGROW>
        if hasEdgeAttr, loadVars{end+1} = 'EA_cell'; end %#ok<AGROW>
        if hasGlobals
            loadVars = [loadVars, {'TR_vec', 'ANG_vec'}]; %#ok<AGROW>
        end
        if hasVerLabel
            loadVars{end+1} = 'schema_version_label'; %#ok<AGROW>
        end
        agg = load(aggFile, loadVars{:});
        aggMap = containers.Map(agg.sample_ids, num2cell(1:numel(agg.sample_ids)));
        if hasVerLabel
            versionLabel = char(agg.schema_version_label);
        elseif hasEdgeAttr
            versionLabel = 'V2';
        else
            versionLabel = 'V1';
        end
    end

    X_cell  = cell(G, 1);
    ei_cell = cell(G, 1);
    EA_cell = cell(G, 1);
    N_vec   = zeros(G, 1);
    G_mat   = zeros(G, 2);
    Dense   = [];
    denseSeen = false(G, 1);

    for g = 1:G
        sid = sample_ids{g};
        if ~isKey(aggMap, sid)
            error('load_hybrid_graph_dataset:missing', ...
                'sample_id "%s" not found in graphs_all.mat. Re-run step9.', sid);
        end
        idx = aggMap(sid);
        X_cell{g} = ensure_boundary_feature(single(agg.X_cell{idx}));
        ei_cell{g} = double(agg.ei_cell{idx});
        N_vec(g) = agg.N_vec(idx);
        if hasEdgeAttr
            EA_cell{g} = single(agg.EA_cell{idx});
        else
            EA_cell{g} = zeros(0, 0, 'single');
        end
        if hasGlobals
            G_mat(g, :) = [agg.TR_vec(idx), agg.ANG_vec(idx)];
        else
            G_mat(g, :) = parse_sample_id(sid);
        end
    end

    if hasDense
        firstDense = agg.Dense_cell{aggMap(sample_ids{1})};
        Dense = zeros(size(firstDense, 1), size(firstDense, 2), size(firstDense, 3), G, 'uint8');
        for g = 1:G
            Dense(:, :, :, g) = uint8(agg.Dense_cell{aggMap(sample_ids{g})});
        end
    else
        Dense = rasterize_from_graphs(X_cell, 128);
    end

    datasetVersion = versionLabel;
else
    fprintf('load_hybrid_graph_dataset: graphs_all.mat not found - loading %d files individually.\n', G);
    X_cell  = cell(G, 1);
    ei_cell = cell(G, 1);
    EA_cell = cell(G, 1);
    N_vec   = zeros(G, 1);
    G_mat   = zeros(G, 2);
    Dense   = [];
    sawEdgeAttr   = false;
    sawVerLabel   = '';

    for g = 1:G
        sid = sample_ids{g};
        fpath = fullfile(samplesDir, sid, 'sample.mat');
        if ~isfile(fpath)
            error('load_hybrid_graph_dataset:missing', 'sample.mat not found: %s', fpath);
        end
        d = load(fpath);
        gd = d.gnn_data;
        X_cell{g} = ensure_boundary_feature(single(gd.x.'));
        ei_cell{g} = double(gd.edge_index);
        N_vec(g) = gd.num_nodes;
        G_mat(g, :) = [gd.tr_ratio, gd.ang_deg];
        if isfield(gd, 'edge_attr') && ~isempty(gd.edge_attr)
            EA_cell{g} = single(gd.edge_attr.');
            sawEdgeAttr = true;
        else
            EA_cell{g} = zeros(0, 0, 'single');
        end
        if isfield(gd, 'schema_version_label') && isempty(sawVerLabel)
            sawVerLabel = char(gd.schema_version_label);
        end
        if isfield(d, 'dense_data') && isfield(d.dense_data, 'raster')
            if isempty(Dense)
                r = d.dense_data.raster;
                Dense = zeros(size(r, 1), size(r, 2), size(r, 3), G, 'uint8');
            end
            Dense(:, :, :, g) = uint8(d.dense_data.raster);
            denseSeen(g) = true;
        end
    end

    if isempty(Dense)
        Dense = rasterize_from_graphs(X_cell, 128);
    elseif any(~denseSeen)
        fallbackDense = rasterize_from_graphs(X_cell, size(Dense, 1));
        Dense(:, :, :, ~denseSeen) = fallbackDense(:, :, :, ~denseSeen);
    end

    if ~isempty(sawVerLabel)
        datasetVersion = sawVerLabel;
    elseif sawEdgeAttr
        datasetVersion = 'V2';
    else
        datasetVersion = 'V1';
    end
end

N_vec = double(N_vec);
end

function X = ensure_boundary_feature(X)
if size(X, 1) < 4
    X = [X; zeros(1, size(X, 2), 'single')];
end
end

function Dense = rasterize_from_graphs(X_cell, gridSize)
G = numel(X_cell);
Dense = zeros(gridSize, gridSize, 2, G, 'uint8');
for g = 1:G
    X = X_cell{g};
    xy = double(X(1:2, :).');
    if size(X, 1) >= 7
        % V2 one-hot: col 4 is interior; any of cols 5..7 is "boundary".
        boundary = any(X(5:7, :), 1).';
    else
        boundary = logical(X(4, :).');
    end
    xLo = min(xy(:, 1)); xHi = max(xy(:, 1));
    yLo = min(xy(:, 2)); yHi = max(xy(:, 2));
    cols = 1 + round((xy(:, 1) - xLo) / max(xHi - xLo, eps) * (gridSize - 1));
    rows = 1 + round((xy(:, 2) - yLo) / max(yHi - yLo, eps) * (gridSize - 1));
    cols = min(gridSize, max(1, cols));
    rows = min(gridSize, max(1, rows));
    idx = sub2ind([gridSize, gridSize], rows, cols);
    Dense(idx + (g - 1) * gridSize * gridSize * 2) = 1;
    bidx = idx(boundary);
    Dense(bidx + gridSize * gridSize + (g - 1) * gridSize * gridSize * 2) = 1;
end
end

function out = parse_sample_id(sample_id)
tr_ratio = NaN;
ang_deg = NaN;
tk = regexp(sample_id, 'tr(\d+)', 'tokens', 'once');
if ~isempty(tk), tr_ratio = str2double(tk{1}) / 100; end
tk = regexp(sample_id, 'ang(\d+)', 'tokens', 'once');
if ~isempty(tk), ang_deg = str2double(tk{1}); end
out = [tr_ratio, ang_deg];
end
