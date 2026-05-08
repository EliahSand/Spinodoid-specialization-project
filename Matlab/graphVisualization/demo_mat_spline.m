%% Demo: visualize one MAT/spline v4 GNN graph sample
% Produces the same graphVisualization tabbed view as demo.m, but reads
% schema v4 samples from dataset_mat_spline.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
gnnRoot = fullfile(fileparts(scriptDir), 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset_mat_spline', 'samples');

if ~isfolder(samplesDir)
    error('MAT/spline samples directory not found: %s', samplesDir);
end

sample_id = first_sample_id(samplesDir);
visualize_mat_spline_graph(sample_id, samplesDir);

function sample_id = first_sample_id(samplesDir)
items = dir(samplesDir);
for k = 1:numel(items)
    if items(k).isdir && items(k).name(1) ~= '.' && ...
            isfile(fullfile(samplesDir, items(k).name, 'sample.mat'))
        sample_id = items(k).name;
        return;
    end
end
error('No MAT/spline sample.mat files found under: %s', samplesDir);
end
