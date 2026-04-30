function [X_pad, nodeMask, normStats] = pad_and_normalize_hybrid_graphs(X_cell, N_vec, trainMask, normStats)
%PAD_AND_NORMALIZE_HYBRID_GRAPHS Pad structural graphs for the hybrid model.
%   Node features are [x, y, radius, boundary]. Coordinates are normalized
%   per graph to preserve shape in a common [-1,1] frame; radius is z-scored
%   with train-split statistics only.

G = numel(N_vec);
maxN = max(N_vec);
F = 4;

if nargin < 4 || isempty(normStats)
    rVals = [];
    for g = 1:G
        if trainMask(g)
            rVals = [rVals, X_cell{g}(3, :)]; %#ok<AGROW>
        end
    end
    normStats.r_mean = mean(rVals);
    normStats.r_std = max(std(rVals), 1e-8);
end

X_pad = zeros(F, maxN, G, 'single');
nodeMask = false(1, maxN, G);

for g = 1:G
    x = X_cell{g};
    if size(x, 1) < 4
        x = [x; zeros(1, size(x, 2), 'single')];
    end
    Ng = N_vec(g);

    for dim = 1:2
        lo = min(x(dim, :));
        hi = max(x(dim, :));
        x(dim, :) = 2 * ((x(dim, :) - lo) / max(hi - lo, eps)) - 1;
    end
    x(3, :) = (x(3, :) - normStats.r_mean) / normStats.r_std;
    x(4, :) = single(x(4, :) > 0);

    X_pad(:, 1:Ng, g) = x(1:F, :);
    nodeMask(1, 1:Ng, g) = true;
end
end
