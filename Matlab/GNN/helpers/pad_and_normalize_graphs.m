function [X_pad, nodeMask, normStats] = pad_and_normalize_graphs(X_cell, N_vec, trainMask, normStats)
%PAD_AND_NORMALIZE_GRAPHS Pad to maxN and normalize node features.
%   Global z-score for x,y,radius using train-split statistics.
%   trainMask  — logical G×1, true for training samples.
%   normStats  — (optional) precomputed stats; if omitted, computed from train split.

G    = numel(N_vec);
maxN = max(N_vec);
F    = size(X_cell{1}, 1);   % x, y, radius, boundary

% Coordinate and radius stats from training graphs (computed before padding)
if nargin < 4 || isempty(normStats)
    x_vals = []; y_vals = []; r_vals = [];
    for g = 1:G
        if trainMask(g)
            x_vals = [x_vals, X_cell{g}(1,:)]; %#ok<AGROW>
            y_vals = [y_vals, X_cell{g}(2,:)]; %#ok<AGROW>
            r_vals = [r_vals, X_cell{g}(3,:)]; %#ok<AGROW>
        end
    end
    normStats.x_mean = mean(x_vals);
    normStats.x_std  = max(std(x_vals), 1e-8);
    normStats.y_mean = mean(y_vals);
    normStats.y_std  = max(std(y_vals), 1e-8);
    normStats.r_mean = mean(r_vals);
    normStats.r_std  = max(std(r_vals), 1e-8);
end

X_pad    = zeros(F, maxN, G, 'single');
nodeMask = false(1, maxN, G);

for g = 1:G
    x  = X_cell{g};   % F×Ng
    Ng = N_vec(g);

    % Global z-score x,y using train stats (domain is fixed across all samples)
    x(1,:) = (x(1,:) - normStats.x_mean) / normStats.x_std;
    x(2,:) = (x(2,:) - normStats.y_mean) / normStats.y_std;

    % Z-score radius using train stats
    x(3,:) = (x(3,:) - normStats.r_mean) / normStats.r_std;

    X_pad(:, 1:Ng, g)    = x;
    nodeMask(1, 1:Ng, g) = true;
end
end