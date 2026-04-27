function Zhat = gcn_forward(params, X, A_hat, K, nodeMask, gpuUse)
%GCN_FORWARD Graph-level GCN: K residual message-passing layers + masked mean-pool.
%   X:        F×maxN×G (plain or dlarray)
%   A_hat:    G×1 cell of sparse maxN×maxN normalized adjacency
%   nodeMask: 1×maxN×G logical or single
%   Zhat:     nComp×1×G dlarray

B = size(X, 3);

% Pre-convert adjacency to dense single (optionally GPU)
A_dense = cell(B, 1);
for b = 1:B
    A = single(full(A_hat{b}));
    if gpuUse, A = gpuArray(A); end
    A_dense{b} = A;
end

% Embedding layer
u = tanh(pagemtimes(params.Embedding.W, X) + params.Embedding.b);

% GCN layers with residual connection
for i = 1:K
    v    = pagemtimes(params.("GNN_"+i).W, u) + params.("GNN_"+i).b;
    Vagg = zeros(size(v), 'like', v);
    for b = 1:B
        Vagg(:,:,b) = dlarray((A_dense{b} * v(:,:,b).').');
    end
    u = u + tanh(Vagg);
end

% Masked mean-pool: F×maxN×G → F×1×G
if isa(u, 'dlarray')
    uData = extractdata(u);
else
    uData = u;
end
mask = cast(nodeMask, 'like', uData);
mask = dlarray(mask);           % 1×maxN×G, no gradient

u_m = u .* mask;
h   = sum(u_m, 2) ./ sum(mask, 2);   % F×1×G (sum(mask,2) > 0 for all valid graphs)

% Pool MLP → Decoder
h    = tanh(pagemtimes(params.PoolMLP.W, h) + params.PoolMLP.b);
Zhat = pagemtimes(params.Decoder.W, h) + params.Decoder.b;   % nComp×1×G
end
