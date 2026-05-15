function [X_pad, nodeMask, normStats] = pad_and_normalize_hybrid_graphs(X_cell, N_vec, trainMask, normStats)
%PAD_AND_NORMALIZE_HYBRID_GRAPHS Pad structural graphs for the hybrid model.
%   V1 layout: F=4 — [x, y, radius, boundary_flag]. Boundary thresholded to {0,1}.
%   V2 layout: F=7 — [x, y, radius, bc_oneHot4]. One-hot already in {0,1}; pass through.
%   Coordinates are box-normalized per graph; radius is z-scored with train stats.

G    = numel(N_vec);
maxN = max(N_vec);
F    = size(X_cell{1}, 1);

if nargin < 4 || isempty(normStats)
    rVals = [];
    for g = 1:G
        if trainMask(g)
            rVals = [rVals, X_cell{g}(3, :)]; %#ok<AGROW>
        end
    end
    normStats.r_mean = mean(rVals);
    normStats.r_std  = max(std(rVals), 1e-8);
end

X_pad    = zeros(F, maxN, G, 'single');
nodeMask = false(1, maxN, G);

for g = 1:G
    x  = X_cell{g};
    Ng = N_vec(g);

    if size(x, 1) ~= F
        error('pad_and_normalize_hybrid_graphs:inconsistentF', ...
            'Graph %d has %d features but expected %d.', g, size(x, 1), F);
    end

    for dim = 1:2
        lo = min(x(dim, :));
        hi = max(x(dim, :));
        x(dim, :) = (x(dim, :) - lo) / max(hi - lo, eps);
    end
    x(3, :) = (x(3, :) - normStats.r_mean) / normStats.r_std;

    if F == 4
        x(4, :) = single(x(4, :) > 0);
    end
    % For V2 (F == 7), cols 4..7 are already a {0,1} one-hot; pass through.

    X_pad(:, 1:Ng, g)    = x(1:F, :);
    nodeMask(1, 1:Ng, g) = true;
end
end
