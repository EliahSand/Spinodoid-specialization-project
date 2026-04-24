function [loss, grad] = model_loss(params, X, A_hat, K, nodeMask, gpuUse, Z_target, w, G_global)
%MODEL_LOSS Weighted mean squared error in standardized PCA space.
%   w: nComp×1×1 per-component weight (default = uniform). Use w = Z_std.^2
%   to recover raw-Z MSE (i.e. u3-space MSE up to the orthonormal PCA basis).
%   G_global: nGlobal×1×B dlarray of standardized global features (optional).

if nargin < 8 || isempty(w), w = 1; end
if nargin < 9, G_global = []; end

Zhat = gcn_forward(params, X, A_hat, K, nodeMask, gpuUse, G_global);
err  = (Zhat - Z_target).^2;              % nComp × 1 × B
loss = mean(sum(w .* err, 1), 'all');     % weighted over components, mean over batch

if nargout > 1
    grad = dlgradient(loss, params);
end
end