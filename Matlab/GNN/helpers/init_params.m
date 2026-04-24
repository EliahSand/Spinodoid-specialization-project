function params = init_params(F, hiddenDim, poolDim, nComp, K, nGlobal, readoutDim)
%INIT_PARAMS Build GCN parameter struct with Glorot-initialized weights.
%   Architecture: Embedding → K×GCN → mean-pool → concat(globals) → Readout → Decoder
%   nGlobal: number of global graph features concatenated after the pool (default 0).
%   readoutDim: dimension of readout head after pooling (default = hiddenDim).

if nargin < 6, nGlobal = 0; end
if nargin < 7 || isempty(readoutDim), readoutDim = hiddenDim; end

% Embedding layer
params.Embedding.W = initializeGlorot([hiddenDim, F]);
params.Embedding.b = zeros(hiddenDim, 1);

% GCN layers with residual connections
for i = 1:K
    params.("GNN_"+i).W = initializeGlorot([hiddenDim, hiddenDim]);
    params.("GNN_"+i).b = zeros(hiddenDim, 1);
end

% Readout head (takes GCN output + concatenated globals)
params.Readout1.W = initializeGlorot([readoutDim, hiddenDim + nGlobal]);
params.Readout1.b = zeros(readoutDim, 1);
params.Readout2.W = initializeGlorot([readoutDim, readoutDim]);
params.Readout2.b = zeros(readoutDim, 1);

% Output decoder
params.Decoder.W = initializeGlorot([nComp, readoutDim]);
params.Decoder.b = zeros(nComp, 1);

params = dlupdate(@dlarray, params);
end