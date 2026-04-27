function D = prepare_dense_raster_inputs(Dense)
%PREPARE_DENSE_RASTER_INPUTS Add coordinate channels to uint8 dense rasters.
%   Dense: H x W x 2 x G uint8 [occupancy, boundary]
%   D:     H x W x 4 x G single [occupancy, boundary, x_grid, y_grid]

H = size(Dense, 1);
W = size(Dense, 2);
G = size(Dense, 4);

D = zeros(H, W, 4, G, 'single');
D(:, :, 1:2, :) = single(Dense(:, :, 1:2, :));

[xx, yy] = meshgrid(single(linspace(0, 1, W)), single(linspace(0, 1, H)));
D(:, :, 3, :) = repmat(xx, 1, 1, 1, G);
D(:, :, 4, :) = repmat(yy, 1, 1, 1, G);
end
