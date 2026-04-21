function A_hat = build_norm_adjacency(ei_cell, N_vec, maxN)
%BUILD_NORM_ADJACENCY Compute D^{-1/2}(A+I)D^{-1/2} for each graph, padded to maxN.
%   ei_cell{g} — 2×Eg double, 1-based node indices
%   N_vec(g)   — number of valid nodes
%   maxN       — padded graph size (shared across batch)

G     = numel(N_vec);
A_hat = cell(G, 1);

for g = 1:G
    Ng = N_vec(g);
    ei = ei_cell{g};      % 2×Eg

    if isempty(ei)
        A = sparse(Ng, Ng);
    else
        rows = [ei(1,:), ei(2,:)];
        cols = [ei(2,:), ei(1,:)];
        A    = spones(sparse(rows, cols, ones(1, 2*size(ei,2)), Ng, Ng));
    end

    % Embed in maxN×maxN and add self-loops
    Ap = spalloc(maxN, maxN, nnz(A) + maxN);
    Ap(1:Ng, 1:Ng) = A;
    Ap = Ap + speye(maxN);

    % Symmetric normalization: D^{-1/2} (A+I) D^{-1/2}
    d  = full(sum(Ap, 2));
    d(d == 0) = 1;
    Dm = spdiags(d .^ (-0.5), 0, maxN, maxN);
    A_hat{g} = Dm * Ap * Dm;
end
end
