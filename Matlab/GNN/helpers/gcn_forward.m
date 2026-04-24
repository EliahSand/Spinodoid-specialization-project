function Zhat = gcn_forward(params, X, A_hat, K, nodeMask, gpuUse, G_global)
%GCN_FORWARD Graph-level GCN with mean pooling and concatenated global features.
%   X:        F×maxN×B (plain or dlarray)
%   A_hat:    B×1 cell of sparse maxN×maxN normalized adjacency
%   nodeMask: 1×maxN×B logical or single
%   G_global: nGlobal×1×B dlarray of standardized global features (tr_ratio, ang_deg)
%   Zhat:     nComp×1×B dlarray

B = size(X, 3);

% Pre-convert adjacency to dense single (optionally GPU)
A_dense = cell(B, 1);
for b = 1:B
    A = single(full(A_hat{b}));
    if gpuUse, A = gpuArray(A); end
    A_dense{b} = A;
end

% Prepare mask for GCN layers (must zero out padding activations)
mask = dlarray(single(nodeMask));  % 1×maxN×B, no gradient

% Embedding layer
u = relu(pagemtimes(params.Embedding.W, X) + params.Embedding.b);
u = u .* mask;                  % Zero out padded nodes

% GCN layers with residual connection
for i = 1:K
    v    = pagemtimes(params.("GNN_"+i).W, u) + params.("GNN_"+i).b;
    Vagg = zeros(size(v), 'like', v);
    for b = 1:B
        Vagg(:,:,b) = dlarray((A_dense{b} * v(:,:,b).').');
    end
    u = relu(u + Vagg);
    u = u .* mask;              % Zero out padded nodes after each layer
end

% Masked mean-pool: extract graph-level structure summary
h_pool_gcn = sum(u .* mask, 2) ./ (sum(mask, 2) + eps);  % hiddenDim×1×B

% Concatenate global features (loading context)
if nargin >= 7 && ~isempty(G_global)
    if ~isa(G_global, 'dlarray'), G_global = dlarray(G_global); end
    h_pool = cat(1, h_pool_gcn, G_global);  % (hiddenDim+nGlobal)×1×B
else
    h_pool = h_pool_gcn;
end

% Readout head
h1 = relu(pagemtimes(params.Readout1.W, h_pool) + params.Readout1.b);  % readoutDim×1×B
h2 = relu(pagemtimes(params.Readout2.W, h1) + params.Readout2.b);      % readoutDim×1×B

% Output decoder
Zhat = pagemtimes(params.Decoder.W, h2) + params.Decoder.b;   % nComp×1×B
end