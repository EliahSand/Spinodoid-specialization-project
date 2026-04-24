function visualize_graph(sample_id, samplesDir, targetsDir, options)
%VISUALIZE_GRAPH Visualize structural graph and analyze features
%   sample_id:   Sample ID string (e.g., 'sheetCone_tr100_ang000_lamellar_N128_1x1')
%   samplesDir:  Path to samples directory (e.g., '.../data/dataset/samples')
%   targetsDir:  Path to targets directory (e.g., '.../data/dataset/targets')
%   options:     Optional struct with fields:
%       - showDeformation: false (default) or true (requires U3 data)
%       - figName: figure name for saving
%       - compareMode: 'none' (default), 'angle', or 'ratio'
%       - fixedValue: fixed tr_ratio or ang_deg for comparison
%       - nSamples: number of samples for comparison (default 4)

if nargin < 4 || isempty(options)
    options = struct('showDeformation', false, 'figName', 'graph_analysis', ...
                     'compareMode', 'none', 'fixedValue', 30, 'nSamples', 4);
end

% Add GNN helpers to path
gnnRoot = fileparts(fileparts(samplesDir));
addpath(genpath(fullfile(gnnRoot, 'helpers')));

%% Load graph data

try
    [X_cell, ei_cell, N_vec, G_mat] = load_graph_dataset({sample_id}, samplesDir);
    X = X_cell{1};          % 4 × N (x, y, radius, boundary)
    ei = ei_cell{1};        % 2 × E
    N = N_vec(1);
    tr_ratio = G_mat(1,1);
    ang_deg = G_mat(1,2);
catch ME
    error('Failed to load graph: %s', ME.message);
end

%% Load PCA targets (if available and requested)

pca_target = [];
if options.showDeformation
    try
        t = load(fullfile(targetsDir, 'pca_targets.mat'), 'Z_train', 'train_ids');
        [~, idx] = ismember(sample_id, t.train_ids);
        if idx > 0
            pca_target = t.Z_train(idx, :);
        end
    catch
        warning('Could not load PCA targets');
    end
end

%% Create figure

if strcmp(options.compareMode, 'none')
    h = figure('Name', options.figName, 'Position', [100 100 1200 800]);

    % Plot 1: Graph topology
    subplot(2, 3, 1);
    plot_graph_topology(X, ei, sprintf('Topology: tr=%.2f, ang=%d°', tr_ratio, ang_deg));

    % Plot 2: Node feature distributions
    subplot(2, 3, 2);
    plot_feature_distributions(X);

    % Plot 3: Scatter plot: position vs radius
    subplot(2, 3, 3);
    plot_position_radius(X);

    % Plot 4: Degree distribution
    subplot(2, 3, 4);
    plot_degree_distribution(ei, N);

    % Plot 5: Boundary nodes
    subplot(2, 3, 5);
    plot_boundary_nodes(X, ei, sprintf('Boundary nodes: %.1f%% of nodes', sum(X(4,:))/N*100));

    % Plot 6: Edge length distribution
    subplot(2, 3, 6);
    plot_edge_lengths(X, ei);

    sgtitle(sprintf('Graph Analysis: %s', sample_id));

else
    % Comparison mode
    visualize_comparison(samplesDir, targetsDir, options);
end

end

%% Helper functions

function plot_graph_topology(X, ei, title_str)
    x = X(1, :);
    y = X(2, :);
    radius = X(3, :);
    boundary = double(X(4, :));

    % Plot edges
    for e = 1:size(ei, 2)
        n1 = ei(1, e);
        n2 = ei(2, e);
        line([x(n1), x(n2)], [y(n1), y(n2)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end

    % Plot nodes
    hold on;
    if (max(radius) - min(radius)) > eps
        norm_radius = (radius - min(radius)) / (max(radius) - min(radius));
        node_size = 20 + norm_radius * 80;  % 20-100
    else
        node_size = 50;  % All same size if radius is constant
    end

    scatter(x, y, node_size, boundary + 1, 'filled', 'MarkerFaceAlpha', 0.7);
    colormap(jet);
    c = colorbar;
    c.Ticks = [1 2];
    c.TickLabels = {'Interior', 'Boundary'};
    c.Label.String = 'Boundary';
    xlabel('x'); ylabel('y');
    title(title_str);
    axis equal tight;
    hold off;
end

function plot_feature_distributions(X)
    % Radius distribution
    subplot(2, 2, 1);
    histogram(X(3, :), 30);
    xlabel('Radius'); ylabel('Count');
    title('Radius Distribution');

    % Boundary count
    subplot(2, 2, 2);
    pie([sum(X(4,:)), numel(X(4,:))-sum(X(4,:))], {'Boundary', 'Interior'});
    title('Node Types');

    % X coordinate
    subplot(2, 2, 3);
    histogram(X(1, :), 30);
    xlabel('x position'); ylabel('Count');

    % Y coordinate
    subplot(2, 2, 4);
    histogram(X(2, :), 30);
    xlabel('y position'); ylabel('Count');
end

function plot_position_radius(X)
    scatter(X(1, :), X(2, :), 20, X(3, :), 'filled', 'MarkerFaceAlpha', 0.7);
    xlabel('x'); ylabel('y');
    colorbar;
    title('Position colored by Radius');
    colormap(jet);
end

function plot_degree_distribution(ei, N)
    deg = zeros(N, 1);
    for e = 1:size(ei, 2)
        deg(ei(1, e)) = deg(ei(1, e)) + 1;
        deg(ei(2, e)) = deg(ei(2, e)) + 1;
    end

    histogram(deg, 0:max(deg)+1);
    xlabel('Degree'); ylabel('Count');
    title(sprintf('Degree Distribution (mean=%.2f)', mean(deg)));
end

function plot_boundary_nodes(X, ei, title_str)
    x = X(1, :);
    y = X(2, :);
    boundary = X(4, :);

    % Plot edges
    for e = 1:size(ei, 2)
        n1 = ei(1, e);
        n2 = ei(2, e);
        line([x(n1), x(n2)], [x(n1), x(n2)], [y(n1), y(n2)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end

    % Plot nodes (color by boundary)
    hold on;
    scatter(x(boundary==0), y(boundary==0), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(x(boundary==1), y(boundary==1), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
    legend('Interior', 'Boundary', 'Location', 'best');
    xlabel('x'); ylabel('y');
    title(title_str);
    axis equal tight;
    hold off;
end

function plot_edge_lengths(X, ei)
    lengths = [];
    for e = 1:size(ei, 2)
        n1 = ei(1, e);
        n2 = ei(2, e);
        dx = X(1, n1) - X(1, n2);
        dy = X(2, n1) - X(2, n2);
        lengths = [lengths; sqrt(dx^2 + dy^2)];
    end

    histogram(lengths, 30);
    xlabel('Edge Length'); ylabel('Count');
    title(sprintf('Edge Length Distribution (mean=%.4f)', mean(lengths)));
end

function visualize_comparison(samplesDir, targetsDir, options)
    % Implementation for comparison mode (placeholder)
    % Would load multiple samples and compare them

    % Get all sample IDs
    files = dir(fullfile(samplesDir, '*.mat'));
    sample_ids = {files.name};
    sample_ids = strrep(sample_ids, '.mat', '');

    % Filter by fixed value
    fixed = options.fixedValue;
    if strcmp(options.compareMode, 'angle')
        mask = contains(sample_ids, sprintf('_ang%03d_', round(fix(fixed))));
    else
        mask = contains(sample_ids, sprintf('_tr%03d_', round(fix(fixed))));
    end
    filtered = sample_ids(mask);
    filtered = filtered(1:min(options.nSamples, numel(filtered)));

    % Plot comparison
    n = numel(filtered);
    h = figure('Name', options.figName, 'Position', [100 100 800*n, 800]);

    for i = 1:n
        subplot(1, n, i);
        try
            [X_cell, ei_cell, ~, G_mat] = load_graph_dataset({filtered{i}}, samplesDir);
            X = X_cell{1};
            ei = ei_cell{1};
            tr = G_mat(1,1);
            ang = G_mat(1,2);

            plot_graph_topology(X, ei, sprintf('tr=%.2f, ang=%d°', tr, ang));
        catch
            title('Failed to load');
        end
    end

    if strcmp(options.compareMode, 'angle')
        title(sprintf('Varying Angle (tr=%s)', options.fixedValue));
    else
        title(sprintf('Varying Thickness Ratio (ang=%s)', options.fixedValue));
    end
end