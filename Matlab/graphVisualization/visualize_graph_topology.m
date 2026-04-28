function [figModel, figFull] = visualize_graph_topology(sample_id, samplesDir, varargin)
%VISUALIZE_GRAPH_TOPOLOGY Compatibility wrapper.
%   The graph visualization now has one purpose: show the GNN input graph
%   and the corresponding full reference graph. Use visualize_graph directly
%   in new code.

if nargin < 2
    samplesDir = [];
end

[figModel, figFull] = visualize_graph(sample_id, samplesDir);
end
