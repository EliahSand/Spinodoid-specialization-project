function plot_overlay_series(overlayData, outDir)
%PLOT_OVERLAY_SERIES Create overlay scatter plots for U1, U2, U3 per thickness ratio.

arguments
    overlayData struct
    outDir (1, :) char
end

ensure_dir(outDir);
fields = fieldnames(overlayData);
colorMap = lines(10);

for fIdx = 1:numel(fields)
    fieldName = fields{fIdx};
    ratios = fieldnames(overlayData.(fieldName));
    for r = 1:numel(ratios)
        trLabel = ratios{r};
        data = overlayData.(fieldName).(trLabel);
        if isempty(data.thetaDeg)
            continue;
        end
        fieldLabel = format_field_label(fieldName);
        ratioLabel = format_ratio_label(trLabel);

        fig = figure('Visible', 'off');
        ax = axes(fig);
        hold(ax, 'on');
        for k = 1:numel(data.thetaDeg)
            solidVals = data.solid{k};
            shellVals = data.shell{k};
            cIdx = mod(k-1, size(colorMap, 1)) + 1;
            scatter(ax, solidVals, shellVals, 18, colorMap(cIdx, :), 'filled', ...
                'DisplayName', sprintf('\\theta = %d°', data.thetaDeg(k)), ...
                'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', colorMap(cIdx, :)*0.6, ...
                'MarkerEdgeAlpha', 0.9);
        end
        allVals = cell2mat([data.solid, data.shell].');
        lims = [min(allVals, [], 'omitnan'), max(allVals, [], 'omitnan')];
        if diff(lims) == 0
            lims = lims + [-1, 1];
        else
            pad = 0.05 * diff(lims);
            lims = lims + [-pad, pad];
        end
        plot(ax, lims, lims, 'k--', 'LineWidth', 1.1, 'DisplayName', 'Ideal agreement');
        xlim(ax, lims);
        ylim(ax, lims);
        grid(ax, 'on');
        axis(ax, 'square');
        xlabel(ax, sprintf('%s (solid)', fieldLabel), 'Interpreter', 'tex');
        ylabel(ax, sprintf('%s (shell)', fieldLabel), 'Interpreter', 'tex');
        title(ax, sprintf('%s overlay — %s', fieldLabel, ratioLabel), ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        legend(ax, 'Location', 'bestoutside');

        exportgraphics(fig, fullfile(outDir, sprintf('overlay_%s_%s.png', lower(fieldName), trLabel)), 'Resolution', 300);
        close(fig);
    end
end
