function params = init_params(F, hiddenDim, poolDim, nComp, K)
%INIT_PARAMS Build GCN parameter struct with Glorot-initialized weights.
%   Architecture: Embedding → K×GCN → mean-pool → PoolMLP → Decoder

params.Embedding.W = initializeGlorot([hiddenDim, F]);
params.Embedding.b = zeros(hiddenDim, 1);

for i = 1:K
    params.("GNN_"+i).W = initializeGlorot([hiddenDim, hiddenDim]);
    params.("GNN_"+i).b = zeros(hiddenDim, 1);
end

params.PoolMLP.W = initializeGlorot([poolDim, hiddenDim]);
params.PoolMLP.b = zeros(poolDim, 1);

params.Decoder.W = initializeGlorot([nComp, poolDim]);
params.Decoder.b = zeros(nComp, 1);

params = dlupdate(@dlarray, params);
end
