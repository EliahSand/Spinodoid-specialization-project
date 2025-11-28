function plot_curvature_profiles(solid, shell, metrics, trLabel, thetaDeg, outDir)
%PLOT_CURVATURE_PROFILES Overlay solid vs shell profiles along the mid-line.
%   Two rows: top = solid vs shell overlay; bottom = absolute error curve.

arguments
    solid table
    shell table
    metrics struct
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

ensure_dir(outDir);

fieldsList = {'U1','U2','U3','S11','S22','SMises'};
[coords, coordLabel, sortIdx] = pick_axis_and_sort(solid);
solid = solid(sortIdx, :);
shell = shell(sortIdx, :);

nFields = numel(fieldsList);
fig = figure('Visible', 'off');
tl = tiledlayout(fig, nFields, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Solid vs shell profiles — %s | \\theta = %d°', trLabel, thetaDeg), 'FontWeight', 'bold');
subtitle(tl, 'Left: overlay; Right: |shell - solid|', 'FontSize', 10);

for i = 1:nFields
    fname = fieldsList{i};
    sVal = double(solid.(fname));
    shVal = double(shell.(fname));
    err = abs(shVal - sVal);

    % Overlay
    ax1 = nexttile(tl);
    plot(ax1, coords, sVal, '-k', 'LineWidth', 1.5, 'DisplayName', 'Solid');
    plot(ax1, coords, shVal, '-r', 'LineWidth', 1.5, 'DisplayName', 'Shell');
    grid(ax1, 'on');
    xlabel(ax1, coordLabel);
    ylabel(ax1, fname);
    m = metrics.(fname);
    title(ax1, sprintf('%s (R=%.3f, MAE=%.3g)', fname, m.R, m.MAE), 'FontSize', 11);
    legend(ax1, 'Location', 'best');

    % Absolute error
    ax2 = nexttile(tl);
    plot(ax2, coords, err, '-b', 'LineWidth', 1.2);
    grid(ax2, 'on');
    xlabel(ax2, coordLabel);
    ylabel(ax2, sprintf('|d%s|', fname));
    title(ax2, sprintf('%s error (Max=%.3g)', fname, max(err, [], 'omitnan')), 'FontSize', 11);
end

exportgraphics(fig, fullfile(outDir, 'profiles_overlay.png'), 'Resolution', 300);
close(fig);

end

function [coord, label, sortIdx] = pick_axis_and_sort(tbl)
coordCandidates = {'X','Y','Z'};
ranges = zeros(1, numel(coordCandidates));
for i = 1:numel(coordCandidates)
    c = coordCandidates{i};
    if ismember(c, tbl.Properties.VariableNames)
        vals = double(tbl.(c));
        ranges(i) = max(vals) - min(vals);
    end
end
[maxRange, idx] = max(ranges);
if maxRange <= 0
    coord = (1:height(tbl))';
    label = 'Node index';
else
    coord = double(tbl.(coordCandidates{idx}));
    label = coordCandidates{idx};
end
[coord, sortIdx] = sort(coord);
end
