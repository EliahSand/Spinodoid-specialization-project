function [loss, grad] = hybrid_model_loss(params, X, A_hat, K, nodeMask, Dense, G_global, gpuUse, Z_target, w, dropoutRate, useGraph, useDense)
%HYBRID_MODEL_LOSS Weighted MSE in standardized PCA space.

if nargin < 10 || isempty(w), w = 1; end
if nargin < 11 || isempty(dropoutRate), dropoutRate = 0; end
if nargin < 12 || isempty(useGraph), useGraph = true; end
if nargin < 13 || isempty(useDense), useDense = true; end

Zhat = hybrid_forward(params, X, A_hat, K, nodeMask, Dense, G_global, gpuUse, dropoutRate, useGraph, useDense);
err = (Zhat - Z_target).^2;
loss = mean(sum(w .* err, 1), 'all');

if nargout > 1
    grad = dlgradient(loss, params);
end
end
