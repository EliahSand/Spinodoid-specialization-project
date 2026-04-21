function [X_pad, nodeMask, normStats] = pad_and_normalize_graphs(X_cell, N_vec, trainMask, normStats)
%PAD_AND_NORMALIZE_GRAPHS Pad to maxN and normalize node features.
%   Per-graph: box-normalize (x,y) to [0,1].
%   Global:    z-score radius using train-split statistics.
%   trainMask  — logical G×1, true for training samples.
%   normStats  — (optional) precomputed stats; if omitted, computed from train split.

G    = numel(N_vec);
maxN = max(N_vec);
F    = size(X_cell{1}, 1);   % 3: x, y, radius

% Radius stats from training graphs (computed before padding)
if nargin < 4 || isempty(normStats)
    r_vals = [];
    for g = 1:G
        if trainMask(g)
            r_vals = [r_vals, X_cell{g}(3,:)]; %#ok<AGROW>
        end
    end
    normStats.r_mean = mean(r_vals);
    normStats.r_std  = max(std(r_vals), 1e-8);
end

X_pad    = zeros(F, maxN, G, 'single');
nodeMask = false(1, maxN, G);

for g = 1:G
    x  = X_cell{g};   % 3×Ng
    Ng = N_vec(g);

    % Box-normalize x,y per graph to [0,1]
    for dim = 1:2
        lo = min(x(dim,:));
        hi = max(x(dim,:));
        x(dim,:) = (x(dim,:) - lo) / max(hi - lo, eps);
    end

    % Z-score radius using train stats
    x(3,:) = (x(3,:) - normStats.r_mean) / normStats.r_std;

    X_pad(:, 1:Ng, g)    = x;
    nodeMask(1, 1:Ng, g) = true;
end
end
