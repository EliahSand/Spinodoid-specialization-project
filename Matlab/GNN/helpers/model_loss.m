function [loss, grad] = model_loss(params, X, A_hat, K, nodeMask, gpuUse, Z_target)
%MODEL_LOSS Mean squared error in standardized PCA space.

Zhat = gcn_forward(params, X, A_hat, K, nodeMask, gpuUse);
loss = mean((Zhat - Z_target).^2, "all");

if nargout > 1
    grad = dlgradient(loss, params);
end
end
