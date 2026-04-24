% visualize_graph_topology.m
%
% Visualize the structural graph for one or more samples.
% Shows node positions, edge connectivity, node radius (color + marker size),
% and boundary vs interior classification.
%
% Usage:
%   Run as script — edit SAMPLE_ID and optional NEIGHBOR_NODE below.
%   Set SHOW_NEIGHBOR_RING = true to highlight a node and its 1-hop ring.

scriptPath = mfilename('fullpath');
scriptDir  = fileparts(scriptPath);
gnnRoot    = fileparts(scriptDir);
helpersDir = fullfile(gnnRoot, 'helpers');
addpath(helpersDir);

samplesDir = fullfile(gnnRoot, 'data', 'dataset', 'samples');

% --- Config -----------------------------------------------------------
SAMPLE_ID         = 'sheetCone_tr100_ang000_lamellar_N128_1x1';
SHOW_NEIGHBOR_RING = true;   % true → pick a node and highlight its ring
NEIGHBOR_NODE      = [];     % [] → auto-pick centroid-nearest node
% ----------------------------------------------------------------------

[X_cell, ei_cell, ~, G_mat] = load_graph_dataset({SAMPLE_ID}, samplesDir);

X  = double(X_cell{1});   % 4 × N  [x; y; radius; boundary]
EI = ei_cell{1};           % 2 × E  (1-based)

x_coords  = X(1,:).';
y_coords  = X(2,:).';
radii     = X(3,:).';
is_bnd    = logical(X(4,:).');

N = numel(x_coords);
E = size(EI, 2);

fprintf('Sample : %s\n', SAMPLE_ID);
fprintf('Nodes  : %d   Edges: %d\n', N, E);
fprintf('Radius : min=%.4f  mean=%.4f  max=%.4f\n', min(radii), mean(radii), max(radii));
fprintf('Boundary nodes: %d / %d\n', sum(is_bnd), N);

% Auto-pick neighbor node: closest to centroid
if SHOW_NEIGHBOR_RING && isempty(NEIGHBOR_NODE)
    cx = mean(x_coords);
    cy = mean(y_coords);
    d  = hypot(x_coords - cx, y_coords - cy);
    [~, NEIGHBOR_NODE] = min(d);
end

% --- Figure 1: Full graph colored by radius ---------------------------
figure('Name', 'Graph Topology — Radius', 'Color', 'w', 'Position', [100 100 900 750]);

% Draw edges (thin grey)
hold on;
for e = 1:E
    s = EI(1,e); t = EI(2,e);
    plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.4);
end

% Marker sizes proportional to radius (normalised to [4, 18] pt)
r_norm    = (radii - min(radii)) / (max(radii) - min(radii) + eps);
markerSz  = 4 + 14 * r_norm;

% Interior nodes colored by radius
int_idx = find(~is_bnd);
scatter(x_coords(int_idx), y_coords(int_idx), markerSz(int_idx).^2, ...
    radii(int_idx), 'filled', 'MarkerEdgeColor', 'none');

% Boundary nodes in red with square marker
bnd_idx = find(is_bnd);
scatter(x_coords(bnd_idx), y_coords(bnd_idx), markerSz(bnd_idx).^2, ...
    'rs', 'filled', 'MarkerEdgeColor', [0.5 0 0], 'LineWidth', 0.5);

cb = colorbar;
cb.Label.String = 'Node radius';
colormap(turbo);
axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('%s\n%d nodes, %d edges', strrep(SAMPLE_ID,'_','\_'), N, E));
legend({'Edges','Interior (radius)','Boundary'}, 'Location','best', 'FontSize', 8);
hold off;

% --- Figure 2: Neighbor ring highlight --------------------------------
if SHOW_NEIGHBOR_RING
    % Find 1-hop neighbors of NEIGHBOR_NODE
    src = EI(1,:); dst = EI(2,:);
    neighbors = unique([dst(src == NEIGHBOR_NODE), src(dst == NEIGHBOR_NODE)]);
    ring_nodes = neighbors(:);

    figure('Name', 'Graph Topology — Neighbor Ring', 'Color', 'w', 'Position', [1020 100 900 750]);
    hold on;

    % All edges (grey)
    for e = 1:E
        s = EI(1,e); t = EI(2,e);
        plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
    end

    % Edges connecting NEIGHBOR_NODE to its ring (blue)
    for e = 1:E
        s = EI(1,e); t = EI(2,e);
        if s == NEIGHBOR_NODE || t == NEIGHBOR_NODE
            plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
                'b-', 'LineWidth', 1.5);
        end
    end

    % All nodes (light, sized by radius)
    scatter(x_coords, y_coords, markerSz.^2, [0.8 0.8 0.8], 'filled');

    % Ring neighbors (cyan)
    scatter(x_coords(ring_nodes), y_coords(ring_nodes), markerSz(ring_nodes).^2 * 1.5, ...
        'c', 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'LineWidth', 1);

    % Center node (red star)
    scatter(x_coords(NEIGHBOR_NODE), y_coords(NEIGHBOR_NODE), 200, ...
        'rp', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

    % Radius labels for center + ring
    for n = [NEIGHBOR_NODE; ring_nodes]
        text(x_coords(n), y_coords(n), sprintf('  r=%.3f', radii(n)), ...
            'FontSize', 6, 'Color', [0.2 0.2 0.2]);
    end

    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('Node %d and its %d 1-hop neighbors\nradius=%.4f  boundary=%d', ...
        NEIGHBOR_NODE, numel(ring_nodes), radii(NEIGHBOR_NODE), is_bnd(NEIGHBOR_NODE)));
    legend({'All edges','Ring edges','All nodes','Ring neighbors','Center node'}, ...
        'Location','best','FontSize',8);
    hold off;
end
