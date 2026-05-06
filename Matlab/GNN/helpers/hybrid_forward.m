function Zhat = hybrid_forward(params, X, A_hat, K, nodeMask, Dense, G_global, gpuUse, dropoutRate, useGraph, useDense)
%HYBRID_FORWARD Dense-raster CNN + structural-GNN forward pass.
%   X:        F x maxN x B
%   Dense:    H x W x 4 x B
%   G_global: nGlobal x 1 x B

if nargin < 9 || isempty(dropoutRate), dropoutRate = 0; end
if nargin < 10 || isempty(useGraph), useGraph = true; end
if nargin < 11 || isempty(useDense), useDense = true; end

if ~isa(G_global, 'dlarray'), G_global = dlarray(G_global); end
B = size(G_global, 3);

if useGraph
    hGraph = graph_branch(params, X, A_hat, K, nodeMask, gpuUse);
else
    graphDim = 2 * size(params.GraphEmbedding.W, 1) * (K + 1);
    hGraph = 0 * repmat(G_global(1, :, :), graphDim, 1, 1);
end

if useDense
    hDense = dense_branch(params, Dense);
else
    cnnDim = 2 * size(params.CNN3.W, 4);
    hDense = 0 * repmat(G_global(1, :, :), cnnDim, 1, 1);
end

if dropoutRate > 0
    keep = dlarray(rand(size(hGraph), 'like', extractdata(hGraph)) > dropoutRate);
    hGraph = hGraph .* keep / (1 - dropoutRate);
    keep = dlarray(rand(size(hDense), 'like', extractdata(hDense)) > dropoutRate);
    hDense = hDense .* keep / (1 - dropoutRate);
end

h = cat(1, hGraph, hDense, G_global);
h = relu(pagemtimes(params.Fusion1.W, h) + params.Fusion1.b);
if dropoutRate > 0
    keep = dlarray(rand(size(h), 'like', extractdata(h)) > dropoutRate);
    h = h .* keep / (1 - dropoutRate);
end
h = relu(pagemtimes(params.Fusion2.W, h) + params.Fusion2.b);
Zhat = pagemtimes(params.Decoder.W, h) + params.Decoder.b;
end

function h = graph_branch(params, X, A_hat, K, nodeMask, gpuUse)
if ~isa(X, 'dlarray'), X = dlarray(X); end
B = size(X, 3);

A_dense = cell(B, 1);
for b = 1:B
    A = single(full(A_hat{b}));
    if gpuUse, A = gpuArray(A); end
    A_dense{b} = A;
end

mask = dlarray(single(nodeMask));
u = relu(pagemtimes(params.GraphEmbedding.W, X) + params.GraphEmbedding.b);
u = u .* mask;

pooled = cell(K + 1, 1);
pooled{1} = graph_pool(u, mask);

for i = 1:K
    v = pagemtimes(params.("Graph_"+i).W, u) + params.("Graph_"+i).b;
    Vagg = zeros(size(v), 'like', v);
    for b = 1:B
        Vagg(:, :, b) = dlarray((A_dense{b} * v(:, :, b).').');
    end
    u = relu(u + Vagg);
    u = u .* mask;
    pooled{i + 1} = graph_pool(u, mask);
end

h = cat(1, pooled{:});
end

function h = graph_pool(u, mask)
hMean = sum(u .* mask, 2) ./ (sum(mask, 2) + eps);
uMasked = u - (1 - cast(mask, 'like', u)) * 1e9;
hMax = max(uMasked, [], 2);
h = cat(1, hMean, hMax);
end

function h = dense_branch(params, Dense)
if ~isa(Dense, 'dlarray'), Dense = dlarray(Dense, 'SSCB'); end

z = dlconv(Dense, params.CNN1.W, params.CNN1.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);
z = dlconv(z, params.CNN2.W, params.CNN2.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);
z = dlconv(z, params.CNN3.W, params.CNN3.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);

B = size(z, 4);
C = size(z, 3);
hMean = stripdims(mean(mean(z, 1), 2));
hMax = stripdims(max(max(z, [], 1), [], 2));
h = cat(1, reshape(hMean, C, 1, B), reshape(hMax, C, 1, B));
end
