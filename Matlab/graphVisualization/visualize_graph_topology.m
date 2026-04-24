function visualize_graph_topology(sample_id, samplesDir, options)
%VISUALIZE_GRAPH_TOPOLOGY Plot structural graph with radius and neighbor ring.
%   sample_id:   Sample ID string
%   samplesDir:  Path to .../data/dataset/samples
%   options:     Optional struct with fields:
%       - showNeighborRing : true (default) — highlight a node's 1-hop ring
%       - neighborNode     : [] (default) — auto-pick centroid-nearest node

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'showNeighborRing'), options.showNeighborRing = true; end
if ~isfield(options, 'neighborNode'),     options.neighborNode = [];       end

gnnRoot    = fileparts(fileparts(samplesDir));
addpath(genpath(fullfile(gnnRoot, 'helpers')));

[X_cell, ei_cell, ~, G_mat] = load_graph_dataset({sample_id}, samplesDir);

X  = double(X_cell{1});   % 4 × N  [x; y; radius; boundary]
EI = ei_cell{1};           % 2 × E

x_coords = X(1,:).';
y_coords = X(2,:).';
radii    = X(3,:).';
is_bnd   = logical(X(4,:).');

N = numel(x_coords);
E = size(EI, 2);

tr_ratio = G_mat(1,1);
ang_deg  = G_mat(1,2);

fprintf('Sample : %s  (tr=%.2f, ang=%g°)\n', sample_id, tr_ratio, ang_deg);
fprintf('Nodes  : %d   Edges: %d\n', N, E);
fprintf('Radius : min=%.4f  mean=%.4f  max=%.4f\n', min(radii), mean(radii), max(radii));
fprintf('Boundary nodes: %d / %d (%.1f%%)\n', sum(is_bnd), N, 100*sum(is_bnd)/N);

r_norm   = (radii - min(radii)) / (max(radii) - min(radii) + eps);
markerSz = 4 + 14 * r_norm;

% --- Figure 1: Full graph colored by radius ---------------------------
figure('Name', sprintf('Radius Topology: %s', sample_id), ...
    'Color', 'w', 'Position', [100 100 900 750]);
hold on;

for e = 1:E
    s = EI(1,e); t = EI(2,e);
    plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.4);
end

int_idx = find(~is_bnd);
scatter(x_coords(int_idx), y_coords(int_idx), markerSz(int_idx).^2, ...
    radii(int_idx), 'filled', 'MarkerEdgeColor', 'none');

bnd_idx = find(is_bnd);
scatter(x_coords(bnd_idx), y_coords(bnd_idx), markerSz(bnd_idx).^2, ...
    'rs', 'filled', 'MarkerEdgeColor', [0.5 0 0], 'LineWidth', 0.5);

cb = colorbar;
cb.Label.String = 'Node radius';
colormap(turbo);
axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('%s\n%d nodes, %d edges | tr=%.2f, ang=%g°', ...
    strrep(sample_id,'_','\_'), N, E, tr_ratio, ang_deg));
legend({'Edges','Interior (radius)','Boundary'}, 'Location','best', 'FontSize', 8);
hold off;

% --- Figure 2: Neighbor ring ------------------------------------------
if options.showNeighborRing
    node = options.neighborNode;
    if isempty(node)
        cx = mean(x_coords); cy = mean(y_coords);
        [~, node] = min(hypot(x_coords - cx, y_coords - cy));
    end

    src = EI(1,:); dst = EI(2,:);
    neighbors = unique([dst(src == node), src(dst == node)]).';

    figure('Name', sprintf('Neighbor Ring: node %d', node), ...
        'Color', 'w', 'Position', [1020 100 900 750]);
    hold on;

    for e = 1:E
        s = EI(1,e); t = EI(2,e);
        plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
    end
    for e = 1:E
        s = EI(1,e); t = EI(2,e);
        if s == node || t == node
            plot([x_coords(s) x_coords(t)], [y_coords(s) y_coords(t)], ...
                'b-', 'LineWidth', 1.5);
        end
    end

    scatter(x_coords, y_coords, markerSz.^2, [0.82 0.82 0.82], 'filled');
    scatter(x_coords(neighbors), y_coords(neighbors), markerSz(neighbors).^2 * 1.5, ...
        'c', 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'LineWidth', 1);
    scatter(x_coords(node), y_coords(node), 200, ...
        'rp', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

    for n = [node; neighbors]
        text(x_coords(n), y_coords(n), sprintf('  r=%.3f', radii(n)), ...
            'FontSize', 6, 'Color', [0.2 0.2 0.2]);
    end

    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('Node %d and its %d 1-hop neighbors (r=%.4f, bnd=%d)', ...
        node, numel(neighbors), radii(node), is_bnd(node)));
    legend({'All edges','Ring edges','All nodes','Ring neighbors','Center'}, ...
        'Location','best', 'FontSize', 8);
    hold off;
end
end
