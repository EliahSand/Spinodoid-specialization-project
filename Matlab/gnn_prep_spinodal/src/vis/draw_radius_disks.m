function draw_radius_disks(ax, x, y, r)
%DRAW_RADIUS_DISKS Draw circle outlines at (x,y) with given radii.
%
%   draw_radius_disks(ax, x, y, r)
%
%   Uses a single vectorised patch call — no per-node loop overhead.

if isempty(x)
    return;
end

faceColor = [0.20, 0.42, 0.85];
faceAlpha = 0.34;
lineWidth = 0.45;

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
    'LineWidth', lineWidth);
end
