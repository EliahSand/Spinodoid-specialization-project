function draw_radius_disks(ax, x, y, r, faceColor, faceAlpha)
%DRAW_RADIUS_DISKS Draw circle outlines at (x,y) with given radii.
%
%   draw_radius_disks(ax, x, y, r, faceColor, faceAlpha)
%
%   Uses a single vectorised patch call — no per-node loop overhead.
%   faceColor  : 1x3 RGB (default [0.20 0.42 0.85])
%   faceAlpha  : edge alpha scalar in [0,1] (default 0.92)

if isempty(x)
    return;
end
if nargin < 5 || isempty(faceColor)
    faceColor = [0.20, 0.42, 0.85];
end
if nargin < 6 || isempty(faceAlpha)
    faceAlpha = 0.92;
end

nNodes = numel(x);
nSeg = 40;
theta = linspace(0, 2 * pi, nSeg + 1);
theta(end) = [];

vertices = zeros(nNodes * nSeg, 2);
faces = reshape(1:(nNodes * nSeg), nSeg, nNodes).';

for i = 1:nNodes
    rows = (i - 1) * nSeg + (1:nSeg);
    vertices(rows, 1) = x(i) + r(i) * cos(theta(:));
    vertices(rows, 2) = y(i) + r(i) * sin(theta(:));
end

patch(ax, ...
    'Faces',     faces, ...
    'Vertices',  vertices, ...
    'FaceColor', 'none', ...
    'EdgeColor', faceColor, ...
    'EdgeAlpha', faceAlpha, ...
    'LineWidth', 0.8);
end
