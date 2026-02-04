function plot_angle_distributions(overlayData, outDir)
%PLOT_ANGLE_DISTRIBUTIONS Compare solid vs shell distributions per angle.

arguments
    overlayData struct
    outDir (1, :) char
end

ensure_dir(outDir);

fieldConfigs = [ ...
    struct('field', 'U3', 'xlabel', '$U_3$ (m)', 'title', 'Out-of-plane displacement', 'bins', 40); ...
    struct('field', 'S_11', 'xlabel', '$S_{11}$ (Pa)', 'title', 'Normal stress $S_{11}$', 'bins', 40); ...
    struct('field', 'S_Mises', 'xlabel', '$\sigma_\mathrm{vM}$ (Pa)', 'title', 'von Mises stress', 'bins', 40)];

solidColor = [0.11 0.35 0.63];
shellColor = [0.80 0.32 0.21];

for cfgIdx = 1:numel(fieldConfigs)
    cfg = fieldConfigs(cfgIdx);
    if ~isfield(overlayData, cfg.field)
        continue;
    end
    ratios = fieldnames(overlayData.(cfg.field));
    for rIdx = 1:numel(ratios)
        trLabel = ratios{rIdx};
        data = overlayData.(cfg.field).(trLabel);
        if isempty(data) || isempty(data.thetaDeg)
            continue;
        end
        ratioLabel = format_ratio_label(trLabel);
        [angles, order] = sort(data.thetaDeg);
        nAngles = numel(angles);
        nCols = min(3, nAngles);
        nRows = ceil(nAngles / nCols);

        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        fieldLabel = format_field_label(cfg.field);
        title(tl, sprintf('%s distribution — %s', fieldLabel, ratioLabel), ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        subtitle(tl, 'Normalized histograms per lamellar angle (solid vs. shell)', 'FontSize', 10);

        for k = 1:nAngles
            idx = order(k);
            solidVals = double(data.solid{idx});
            shellVals = double(data.shell{idx});
            vals = [solidVals(:); shellVals(:)];
            vals = vals(~isnan(vals));
            ax = nexttile(tl);
            if isempty(vals)
                text(ax, 0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
                axis(ax, 'off');
                continue;
            end
            edges = compute_edges(vals, cfg.bins);
            histogram(ax, solidVals, edges, 'Normalization', 'pdf', ...
                'FaceColor', solidColor, 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                'DisplayName', 'Solid');
            hold(ax, 'on');
            histogram(ax, shellVals, edges, 'Normalization', 'pdf', ...
                'FaceColor', shellColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
                'DisplayName', 'Shell');
            muSolid = mean(solidVals, 'omitnan');
            muShell = mean(shellVals, 'omitnan');
            xline(ax, muSolid, '-', 'Color', solidColor, 'LineWidth', 1);
            xline(ax, muShell, '--', 'Color', shellColor, 'LineWidth', 1);
            grid(ax, 'on');
            xlabel(ax, cfg.xlabel, 'Interpreter', 'latex');
            ylabel(ax, 'PDF');
            title(ax, sprintf('\\theta = %d° | \\Delta\\mu = %.2g', ...
                angles(k), muShell - muSolid), 'Interpreter', 'tex');
            legend(ax, 'Location', 'best', 'Box', 'off');
        end

        exportgraphics(fig, fullfile(outDir, sprintf('distribution_%s_%s.png', ...
            lower(cfg.field), trLabel)), 'Resolution', 300);
        close(fig);
    end
end

end

function edges = compute_edges(vals, nBins)
vals = vals(~isnan(vals));
if isempty(vals)
    edges = linspace(-1, 1, max(10, nBins));
    return;
end
mn = min(vals);
mx = max(vals);
if mn == mx
    span = max(abs(mn), 1);
    mn = mn - span;
    mx = mx + span;
end
edges = linspace(mn, mx, max(15, nBins));
end
