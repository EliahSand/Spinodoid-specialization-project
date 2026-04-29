% analyze_pca_relationships.m
%
% Evaluates correlations between the 8 PCA target components.
% Are the components correlated? If so, which pairs are linked?
%
% Inputs:   pca_targets.mat, pca_model.mat
% Outputs:  experiments/pca_relationships/
%             fig1_correlation_matrix.png
%             fig2_pairwise_scatter.png
%             summary.txt

%% ---- Config ---------------------------------------------------------------
REPO_ROOT   = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
TARGETS_DIR = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'data', 'dataset', 'targets');
OUT_DIR     = fullfile(REPO_ROOT, 'Matlab', 'GNN', 'experiments', 'pca_relationships');

%% ---- Load -----------------------------------------------------------------
fprintf('Loading PCA data...\n');
pm = load(fullfile(TARGETS_DIR, 'pca_model.mat'),  'explained');
pt = load(fullfile(TARGETS_DIR, 'pca_targets.mat'), 'Z_train', 'Z_val', 'Z_test');

Z_all    = [pt.Z_train; pt.Z_val; pt.Z_test];  % N_total x K
explained = pm.explained;                        % K x 1
K = size(Z_all, 2);

fprintf('N=%d samples  K=%d components\n', size(Z_all,1), K);

%% ---- Output directory -----------------------------------------------------
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end

%% ---- Figure 1: Correlation matrix -----------------------------------------
R = corrcoef(Z_all);  % K x K

fig = figure('Visible', 'off', 'Position', [0 0 800 700]);
ax = axes(fig);

cmap = interp1([-1 0 1], [0.18 0.38 0.75; 1 1 1; 0.75 0.18 0.18], linspace(-1,1,256));
imagesc(ax, R);
colormap(ax, cmap);
clim(ax, [-1 1]);
colorbar(ax);
axis(ax, 'square');

xticks(ax, 1:K);
yticks(ax, 1:K);
labels = arrayfun(@(k) sprintf('PC%d', k), 1:K, 'UniformOutput', false);
xticklabels(ax, labels);
yticklabels(ax, labels);
title(ax, 'Pearson correlation — PCA scores (all data)', 'FontSize', 13);

for i = 1:K
    for j = 1:K
        val = R(i,j);
        if abs(val) > 0.6
            clr = [1 1 1];
        else
            clr = [0 0 0];
        end
        text(ax, j, i, sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', clr);
    end
end

exportgraphics(fig, fullfile(OUT_DIR, 'fig1_correlation_matrix.png'), 'Resolution', 150);
close(fig);
fprintf('Saved fig1_correlation_matrix.png\n');

%% ---- Figure 2: Pairwise scatter grid --------------------------------------
fig = figure('Visible', 'off', 'Position', [0 0 1400 1400]);

for i = 1:K
    for j = 1:K
        idx = (i-1)*K + j;
        ax = subplot(K, K, idx);

        if i == j
            % Diagonal: histogram
            histogram(ax, Z_all(:,i), 40, 'FaceColor', [0.4 0.6 0.85], 'EdgeColor', 'none');
            title(ax, sprintf('PC%d\n(%.1f%%)', i, explained(i)), 'FontSize', 8);

        elseif i > j
            % Lower triangle: density-coloured scatter
            [N2, xe, ye] = histcounts2(Z_all(:,j), Z_all(:,i), 25);
            xi = discretize(Z_all(:,j), xe);
            yi = discretize(Z_all(:,i), ye);
            valid = xi > 0 & yi > 0;
            density = N2(sub2ind(size(N2), xi(valid), yi(valid)));
            scatter(ax, Z_all(valid,j), Z_all(valid,i), 2, density, 'filled');
            colormap(ax, 'hot');

        else
            % Upper triangle: blank
            axis(ax, 'off');
            continue;
        end

        % Axis labels on edges only
        if i < K
            set(ax, 'XTickLabel', []);
        end
        if j > 1
            set(ax, 'YTickLabel', []);
        end
        if i == K
            xlabel(ax, sprintf('PC%d', j), 'FontSize', 8);
        end
        if j == 1 && i > 1
            ylabel(ax, sprintf('PC%d', i), 'FontSize', 8);
        end

        set(ax, 'FontSize', 7);
    end
end

sgtitle('Pairwise PC score scatter — density coloured (all data)', 'FontSize', 13);
exportgraphics(fig, fullfile(OUT_DIR, 'fig2_pairwise_scatter.png'), 'Resolution', 110);
close(fig);
fprintf('Saved fig2_pairwise_scatter.png\n');

%% ---- Summary (console + summary.txt) --------------------------------------
% Extract off-diagonal pairs sorted by |R|
pairs = [];
for i = 1:K
    for j = i+1:K
        pairs(end+1, :) = [i, j, R(i,j)]; %#ok<AGROW>
    end
end
[~, ord] = sort(abs(pairs(:,3)), 'descend');
pairs_sorted = pairs(ord, :);

% Print to console
fprintf('\n--- Correlation matrix (all data, N=%d) ---\n', size(Z_all,1));
header = sprintf('      ');
for k = 1:K, header = [header sprintf('   PC%d  ', k)]; end %#ok<AGROW>
fprintf('%s\n', header);
for i = 1:K
    row = sprintf('PC%d  ', i);
    for j = 1:K
        row = [row sprintf('%7.4f  ', R(i,j))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\n--- Top-5 strongest off-diagonal correlations ---\n');
for n = 1:min(5, size(pairs_sorted,1))
    fprintf('  PC%d — PC%d :  r = %+.4f\n', ...
        pairs_sorted(n,1), pairs_sorted(n,2), pairs_sorted(n,3));
end

% Write summary.txt
fid = fopen(fullfile(OUT_DIR, 'summary.txt'), 'w');
fprintf(fid, 'PCA Component Correlation Summary\n');
fprintf(fid, 'N = %d samples,  K = %d components\n\n', size(Z_all,1), K);

fprintf(fid, 'Correlation matrix:\n');
fprintf(fid, '%s\n', header);
for i = 1:K
    row = sprintf('PC%d  ', i);
    for j = 1:K
        row = [row sprintf('%7.4f  ', R(i,j))]; %#ok<AGROW>
    end
    fprintf(fid, '%s\n', row);
end

fprintf(fid, '\nTop-5 strongest off-diagonal correlations:\n');
for n = 1:min(5, size(pairs_sorted,1))
    fprintf(fid, '  PC%d - PC%d :  r = %+.4f\n', ...
        pairs_sorted(n,1), pairs_sorted(n,2), pairs_sorted(n,3));
end
fclose(fid);

%% ---- Done -----------------------------------------------------------------
fprintf('\nResults saved to: %s\n', OUT_DIR);
