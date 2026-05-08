function [X_pad, nodeMask, normStats] = pad_and_normalize_hybrid_graphs(X_cell, N_vec, trainMask, normStats)
%PAD_AND_NORMALIZE_HYBRID_GRAPHS Pad structural graphs for the hybrid model.
%   Schema v3 features are [x, y, radius, boundary].
%   Schema v4 MAT/spline features are
%   [x, y, radius, degree, is_endpoint, is_junction, is_internal_control].
%   Coordinates are normalized per graph; radius and v4 degree are z-scored
%   using train-split statistics only.

G = numel(N_vec);
maxN = max(N_vec);
F = max(cellfun(@(x) size(x, 1), X_cell));
F = max(F, 4);

if nargin < 4 || isempty(normStats)
    rVals = [];
    degreeVals = [];
    for g = 1:G
        if trainMask(g)
            xg = X_cell{g};
            if size(xg, 1) >= 3
                rVals = [rVals, xg(3, :)]; %#ok<AGROW>
            end
            if F > 4 && size(xg, 1) >= 4
                degreeVals = [degreeVals, xg(4, :)]; %#ok<AGROW>
            end
        end
    end
    if isempty(rVals)
        rVals = 0;
    end
    normStats.r_mean = mean(rVals);
    normStats.r_std = max(std(rVals), 1e-8);
    normStats.F = F;
    if F > 4
        if isempty(degreeVals)
            degreeVals = 0;
        end
        normStats.degree_mean = mean(degreeVals);
        normStats.degree_std = max(std(degreeVals), 1e-8);
    end
elseif isfield(normStats, 'F')
    F = normStats.F;
end

X_pad = zeros(F, maxN, G, 'single');
nodeMask = false(1, maxN, G);

for g = 1:G
    x = X_cell{g};
    if size(x, 1) < F
        x = [x; zeros(F - size(x, 1), size(x, 2), 'single')];
    end
    Ng = N_vec(g);

    for dim = 1:2
        lo = min(x(dim, :));
        hi = max(x(dim, :));
        x(dim, :) = (x(dim, :) - lo) / max(hi - lo, eps);
    end
    x(3, :) = (x(3, :) - normStats.r_mean) / normStats.r_std;
    if F > 4
        x(4, :) = (x(4, :) - normStats.degree_mean) / normStats.degree_std;
        if F >= 5
            x(5:F, :) = single(x(5:F, :) > 0);
        end
    else
        x(4, :) = single(x(4, :) > 0);
    end

    X_pad(:, 1:Ng, g) = x(1:F, :);
    nodeMask(1, 1:Ng, g) = true;
end
end
