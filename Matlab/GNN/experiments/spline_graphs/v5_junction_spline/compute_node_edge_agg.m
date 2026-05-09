function [EA_pad, normStats] = compute_node_edge_agg(EdgeAttr_cell, ei_cell, N_vec, maxN, trainMask, normStats)
%COMPUTE_NODE_EDGE_AGG  Per-node mean-pooled edge features, z-scored, padded.
%
%   For each graph g, each node i gets a feature vector = mean of
%   EdgeAttr_cell{g}(:, incident_edges). The result is z-scored using
%   train-split statistics (or from normStats when provided) and returned
%   as a padded Fe x maxN x G single array.
%
%   Inputs:
%     EdgeAttr_cell{g} - Fe x Eg single  (output of load_hybrid_graph_dataset)
%     ei_cell{g}       - 2 x Eg double   (edge indices, 1-based)
%     N_vec            - G x 1 node counts
%     maxN             - scalar, max nodes across dataset
%     trainMask        - G x 1 logical
%     normStats        - struct from previous call (or [] to compute fresh)
%
%   Outputs:
%     EA_pad     - Fe x maxN x G single  (padded, z-scored)
%     normStats  - struct with ea_mean (1 x Fe), ea_std (1 x Fe)

G  = numel(N_vec);
Fe = 0;
for g = 1:G
    if ~isempty(EdgeAttr_cell{g})
        Fe = size(EdgeAttr_cell{g}, 1);
        break;
    end
end

EA_pad = zeros(Fe, maxN, G, 'single');
if Fe == 0
    if isempty(normStats)
        normStats = struct('ea_mean', zeros(1,0), 'ea_std', ones(1,0));
    end
    return;
end

% Build per-node aggregated edge features for each graph
agg_cell = cell(G, 1);
for g = 1:G
    Ng = N_vec(g);
    ea = EdgeAttr_cell{g};   % Fe x Eg
    ei = ei_cell{g};         % 2 x Eg
    agg = zeros(Fe, Ng, 'single');
    cnt = zeros(1, Ng, 'single');
    if ~isempty(ea) && ~isempty(ei)
        for e = 1:size(ei, 2)
            u = ei(1,e); v = ei(2,e);
            if u >= 1 && u <= Ng
                agg(:,u) = agg(:,u) + ea(:,e);
                cnt(u)   = cnt(u) + 1;
            end
            if v >= 1 && v <= Ng && v ~= u
                agg(:,v) = agg(:,v) + ea(:,e);
                cnt(v)   = cnt(v) + 1;
            end
        end
        mask = cnt > 0;
        agg(:, mask) = agg(:, mask) ./ cnt(mask);
    end
    agg_cell{g} = agg;
end

% Compute normStats from training samples
if isempty(normStats)
    allVals = [];
    for g = 1:G
        if trainMask(g)
            allVals = [allVals, agg_cell{g}]; %#ok<AGROW>
        end
    end
    if isempty(allVals)
        normStats.ea_mean = zeros(1, Fe, 'single');
        normStats.ea_std  = ones(1, Fe, 'single');
    else
        normStats.ea_mean = mean(allVals, 2).';
        normStats.ea_std  = max(std(allVals, 0, 2), 1e-8).';
    end
end

mu  = normStats.ea_mean(:);   % Fe x 1
sig = normStats.ea_std(:);    % Fe x 1

for g = 1:G
    Ng  = N_vec(g);
    agg = agg_cell{g};
    agg = (agg - mu) ./ sig;
    EA_pad(:, 1:Ng, g) = single(agg);
end
end
