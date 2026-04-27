function metrics = evaluate_hybrid_metrics(Zhat_std, Z_true, Z_mean, Z_std, coeff, u3_mean, U3_true, explained)
%EVALUATE_HYBRID_METRICS Compute PCA-space and decoded u3-space metrics.

Zhat = squeeze(gather(extractdata(Zhat_std))).' .* Z_std + Z_mean;
errZ = Zhat - Z_true;
ssRes = sum(errZ.^2, 1);
ssTot = sum((Z_true - mean(Z_true, 1)).^2, 1);

metrics.perPC_R2 = 1 - ssRes ./ max(ssTot, eps);
metrics.R2_z_unweighted = 1 - sum(ssRes) / max(sum(ssTot), eps);

nComp = size(Z_true, 2);
w = explained(1:nComp) / 100;
metrics.R2_z_weighted = 1 - sum(w(:).' .* ssRes) / max(sum(w(:).' .* ssTot), eps);

U3_pred = Zhat * coeff' + u3_mean;
errU3 = U3_pred - U3_true;
metrics.u3_mse = mean(errU3(:).^2);
metrics.u3_R2 = 1 - sum(errU3(:).^2) / max(sum((U3_true(:) - mean(U3_true(:))).^2), eps);
metrics.Zhat = Zhat;
metrics.U3_pred = U3_pred;
end
