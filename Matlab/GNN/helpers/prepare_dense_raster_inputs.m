function [D, normStats] = prepare_dense_raster_inputs(Dense, trainMask, normStats)
%PREPARE_DENSE_RASTER_INPUTS Add coordinate channels and normalize rasters.
%   Dense: H x W x 2 x G uint8 [occupancy, boundary]
%   D:     H x W x 4 x G single [occupancy, boundary, x_grid, y_grid]
%
%   The binary raster channels are z-scored using train-split pixels only.
%   Coordinate channels are deterministic and centered in [-1, 1].

H = size(Dense, 1);
W = size(Dense, 2);
G = size(Dense, 4);

if nargin < 2 || isempty(trainMask)
    trainMask = true(G, 1);
end
trainMask = logical(trainMask(:));
if numel(trainMask) ~= G
    error('prepare_dense_raster_inputs:badMask', ...
        'trainMask length (%d) must match number of rasters (%d).', numel(trainMask), G);
end
if nargin < 3 || isempty(normStats)
    normStats = struct();
end

D = zeros(H, W, 4, G, 'single');
D(:, :, 1:2, :) = single(Dense(:, :, 1:2, :));

[xx, yy] = meshgrid(single(linspace(-1, 1, W)), single(linspace(-1, 1, H)));
D(:, :, 3, :) = repmat(xx, 1, 1, 1, G);
D(:, :, 4, :) = repmat(yy, 1, 1, 1, G);

if ~isfield(normStats, 'dense_mean') || ~isfield(normStats, 'dense_std')
    dense_mean = zeros(2, 1, 'single');
    dense_std = ones(2, 1, 'single');
    for c = 1:2
        vals = D(:, :, c, trainMask);
        dense_mean(c) = mean(vals(:));
        dense_std(c) = max(std(vals(:)), single(1e-6));
    end
    normStats.dense_mean = dense_mean;
    normStats.dense_std = dense_std;
    normStats.dense_coord_range = single([-1, 1]);
end

for c = 1:2
    D(:, :, c, :) = (D(:, :, c, :) - normStats.dense_mean(c)) ./ normStats.dense_std(c);
end
end
