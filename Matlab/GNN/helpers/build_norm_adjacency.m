function A_hat = build_norm_adjacency(ei_cell, N_vec, maxN, edge_attr_cell, edgeMLP)
%BUILD_NORM_ADJACENCY Symmetric normalized adjacency D^{-1/2}(A+I)D^{-1/2}.
%
%   A_hat = build_norm_adjacency(ei_cell, N_vec, maxN)
%       Binary edges. Returns a cell of sparse maxN×maxN matrices (V1 path).
%
%   A_hat = build_norm_adjacency(ei_cell, N_vec, maxN, edge_attr_cell, edgeMLP)
%       Per-edge scalar weight w_nm = softplus(MLP_edge(edge_attr_nm))
%       built into A_w = w_nm * (A + I)   (self-loop weight = 1)
%       and symmetrically normalized. Returns a cell of dense
%       maxN×maxN dlarray matrices so gradients flow back into edgeMLP.
%
%   ei_cell{g}        — 2×Eg, 1-based node indices
%   edge_attr_cell{g} — Fe×Eg, one column per stored (undirected) edge
%   edgeMLP            — struct with W1 (H×Fe), b1 (H×1), W2 (1×H), b2 (1×1)

G          = numel(N_vec);
A_hat      = cell(G, 1);
useWeighted = nargin >= 5 && ~isempty(edgeMLP) && ~isempty(edge_attr_cell);

for g = 1:G
    Ng = N_vec(g);
    ei = double(ei_cell{g});

    if useWeighted
        A_hat{g} = build_weighted_dense(ei, edge_attr_cell{g}, Ng, maxN, edgeMLP);
    else
        A_hat{g} = build_binary_sparse(ei, Ng, maxN);
    end
end
end


function A_hat = build_binary_sparse(ei, Ng, maxN)
if isempty(ei)
    A = sparse(Ng, Ng);
else
    rows = [ei(1,:), ei(2,:)];
    cols = [ei(2,:), ei(1,:)];
    A    = spones(sparse(rows, cols, ones(1, 2*size(ei,2)), Ng, Ng));
end

Ap = spalloc(maxN, maxN, nnz(A) + maxN);
Ap(1:Ng, 1:Ng) = A;
Ap = Ap + speye(maxN);

d  = full(sum(Ap, 2));
d(d == 0) = 1;
Dm = spdiags(d .^ (-0.5), 0, maxN, maxN);
A_hat = Dm * Ap * Dm;
end


function A_hat = build_weighted_dense(ei, ea, Ng, maxN, edgeMLP) %#ok<INUSL>
% Edge weights w (1×E) via softplus(MLP(ea)), then dense symmetric A_w with
% self-loops of weight 1, then D^{-1/2} A_w D^{-1/2}. Stays dlarray.
E = size(ei, 2);
if E == 0
    % No edges: A = I, normalization is identity. Tie the dtype/device of
    % the result to the MLP weights so downstream casts behave consistently.
    Z = 0 * edgeMLP.W1(1, 1);   % scalar dlarray with right type+device
    A_hat = zeros(maxN, maxN, 'like', Z);
    diagIdx = sub2ind([maxN, maxN], 1:maxN, 1:maxN);
    A_hat(diagIdx) = A_hat(diagIdx) + 1;
    return;
end

if ~isa(ea, 'dlarray'), ea = dlarray(single(ea)); end
if size(ea, 1) ~= size(edgeMLP.W1, 2)
    error('build_norm_adjacency:badEdgeAttr', ...
        'edge_attr has %d rows but MLP expects %d.', size(ea, 1), size(edgeMLP.W1, 2));
end

h = edgeMLP.W1 * ea + edgeMLP.b1;
h = max(h, 0);
z = edgeMLP.W2 * h + edgeMLP.b2;
% Softplus, numerically stable and dlarray-compatible:
%   log(1+e^z) = max(z,0) + log(1 + e^{-|z|})
w = max(z, 0) + log(1 + exp(-abs(z)));   % 1×E


A_w = zeros(maxN, maxN, 'like', w);
src = ei(1, :);
dst = ei(2, :);
idx1 = sub2ind([maxN, maxN], src, dst);
idx2 = sub2ind([maxN, maxN], dst, src);
A_w(idx1) = w;
A_w(idx2) = w;

% Self-loops via diagonal scatter (avoids eye(..,'like',dlarray) which is
% unsupported on some MATLAB versions).
diagIdx = sub2ind([maxN, maxN], 1:maxN, 1:maxN);
A_w(diagIdx) = A_w(diagIdx) + 1;

d        = sum(A_w, 2);                 % maxN×1
d_invs   = 1 ./ sqrt(d + 1e-12);
% (D^{-1/2} A D^{-1/2})_{ij} = A_ij * d_i^{-1/2} * d_j^{-1/2}
A_hat = (A_w .* d_invs) .* d_invs.';
end
