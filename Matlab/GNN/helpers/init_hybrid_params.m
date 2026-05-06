function params = init_hybrid_params(F, hiddenDim, cnnChannels, fusionDim, nComp, K, nGlobal)
%INIT_HYBRID_PARAMS Parameters for CNN + structural GNN hybrid model.

if numel(cnnChannels) ~= 3
    error('init_hybrid_params:badChannels', 'cnnChannels must contain three channel counts.');
end

params.GraphEmbedding.W = initializeGlorot([hiddenDim, F]);
params.GraphEmbedding.b = zeros(hiddenDim, 1);
for i = 1:K
    params.("Graph_"+i).W = initializeGlorot([hiddenDim, hiddenDim]);
    params.("Graph_"+i).b = zeros(hiddenDim, 1);
end

params.CNN1.W = initializeGlorotConv([5, 5, 4, cnnChannels(1)]);
params.CNN1.b = zeros(1, 1, cnnChannels(1));
params.CNN2.W = initializeGlorotConv([3, 3, cnnChannels(1), cnnChannels(2)]);
params.CNN2.b = zeros(1, 1, cnnChannels(2));
params.CNN3.W = initializeGlorotConv([3, 3, cnnChannels(2), cnnChannels(3)]);
params.CNN3.b = zeros(1, 1, cnnChannels(3));

graphDim = 2 * hiddenDim * (K + 1);
cnnDim = 2 * cnnChannels(3);
fusionIn = graphDim + cnnDim + nGlobal;

params.Fusion1.W = initializeGlorot([fusionDim, fusionIn]);
params.Fusion1.b = zeros(fusionDim, 1);
params.Fusion2.W = initializeGlorot([fusionDim, fusionDim]);
params.Fusion2.b = zeros(fusionDim, 1);
params.Decoder.W = initializeGlorot([nComp, fusionDim]);
params.Decoder.b = zeros(nComp, 1);

params = dlupdate(@dlarray, params);
end

function w = initializeGlorotConv(sz)
fanIn = sz(1) * sz(2) * sz(3);
fanOut = sz(1) * sz(2) * sz(4);
scale = sqrt(6 / (fanIn + fanOut));
w = 2 * scale * (rand(sz) - 0.5);
end
