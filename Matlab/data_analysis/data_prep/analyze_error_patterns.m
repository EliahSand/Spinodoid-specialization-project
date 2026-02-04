function analyze_error_patterns(solid, shell, trLabel, thetaDeg, outDir)
%ANALYZE_ERROR_PATTERNS Run PCA on node-wise error features.

arguments
    solid table
    shell table
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

if ~license('test', 'statistics_toolbox') || exist('pca', 'file') ~= 2
    warning('analyze_error_patterns:NoPCA', ...
        'Statistics Toolbox/PCA not available. Skipping PCA for %s at %d°.', trLabel, thetaDeg);
    return;
end

ensure_dir(outDir);

fields = {'U1','U2','U3','S_11','S_22','S_Mises'};
errMat = zeros(height(solid), numel(fields));
for i = 1:numel(fields)
    fname = fields{i};
    errMat(:, i) = shell.(fname) - solid.(fname);
end

coordsX = solid.X;
coordsY = solid.Z;

valid = all(~isnan(errMat), 2);
if ~all(valid)
    warning('analyze_error_patterns:MissingData', ...
        'Excluding %d rows with NaN values from PCA.', sum(~valid));
end
errMat = errMat(valid, :);
coordsX = coordsX(valid);
coordsY = coordsY(valid);

[coeff, score, latent, ~, explained] = pca(errMat, 'Centered', true);

modelTag = sprintf('%s | \\theta = %d°', trLabel, thetaDeg);

% Scree plot
fig1 = figure('Visible', 'off');
bar(explained, 'FaceColor', [0.2 0.4 0.8]);
hold on;
plot(cumsum(explained), '-o', 'LineWidth', 1.2, 'Color', [0.9 0.2 0.2]);
grid on;
xlabel('Principal component');
ylabel('Explained variance (%)');
title(sprintf('PCA scree plot — %s', modelTag), 'FontWeight', 'bold');
legend('Explained variance', 'Cumulative', 'Location', 'best');
exportgraphics(fig1, fullfile(outDir, 'pca_scree.png'), 'Resolution', 300);
close(fig1);

% Spatial map of the first mode + loadings bar chart
fig2 = figure('Visible', 'off');
tl = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Error PCA — %s', modelTag), 'FontWeight', 'bold');
subtitle(tl, 'Spatial PC1 score and component weights');

ax1 = nexttile(tl);
scatter(ax1, coordsX, coordsY, 18, score(:, 1), 'filled', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.85);
axis(ax1, 'equal');
axis(ax1, 'tight');
grid(ax1, 'on');
colormap(ax1, 'turbo');
colorbar(ax1);
xlabel(ax1, 'X');
ylabel(ax1, 'Z');
title(ax1, 'PC1 score (shell-solid)');

ax2 = nexttile(tl);
bar(ax2, coeff(:, 1));
set(ax2, 'XTick', 1:numel(fields), 'XTickLabel', fields);
grid(ax2, 'on');
xlabel(ax2, 'Component');
ylabel(ax2, 'Loading');
title(ax2, 'PC1 loadings');

exportgraphics(fig2, fullfile(outDir, 'pca_mode1.png'), 'Resolution', 300);
close(fig2);
end
