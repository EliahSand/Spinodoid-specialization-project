function plot_aggregated_metrics(summary, trLabels, thetasDeg, outDir)
%PLOT_AGGREGATED_METRICS Visualize metrics across thickness ratios and angles.

arguments
    summary (:,:) struct
    trLabels (1, :) cell
    thetasDeg (1, :) double
    outDir (1, :) char
end

ensure_dir(outDir);

fieldsKey = {'U1','U2','U3','S_11','S_22','S_Mises'};
metricNames = {'MAE', 'RMSE', 'MeanRel'};

nTr = numel(trLabels);
nTheta = numel(thetasDeg);
trVals = cellfun(@(s) sscanf(s, 'tr%d'), trLabels) / 100;
ratioLabels = arrayfun(@(v) format_ratio_value(v), trVals, 'UniformOutput', false);
angleLabels = arrayfun(@(a) sprintf('%d%c', a, char(176)), thetasDeg, 'UniformOutput', false);

data = struct();
for f = 1:numel(fieldsKey)
    field = fieldsKey{f};
    for m = 1:numel(metricNames)
        metric = metricNames{m};
        mat = nan(nTr, nTheta);
        for iTr = 1:nTr
            for jTh = 1:nTheta
                if iTr <= size(summary, 1) && jTh <= size(summary, 2) ...
                        && isfield(summary(iTr, jTh), 'metrics') ...
                        && ~isempty(summary(iTr, jTh).metrics)
                    mat(iTr, jTh) = summary(iTr, jTh).metrics.(field).(metric);
                end
            end
        end
        data.(field).(metric) = mat;
    end
end

for f = 1:numel(fieldsKey)
    field = fieldsKey{f};
    makeLines(field, true, sprintf('agg_%s_vs_ratio.png', lower(field)));
    makeLines(field, false, sprintf('agg_%s_vs_angle.png', lower(field)));
end

heatmapFields = {'U3', 'S_Mises'};
heatmapMetrics = {'MAE', 'MeanRel'};
for hf = 1:numel(heatmapFields)
    field = heatmapFields{hf};
    fig = figure('Visible', 'off');
    tl = tiledlayout(fig, 1, numel(heatmapMetrics), 'TileSpacing', 'compact', 'Padding', 'compact');
    fieldLabel = format_field_label(field);
    title(tl, sprintf('Heatmaps for %s', fieldLabel), 'FontWeight', 'bold', 'Interpreter', 'tex');
    subtitle(tl, 'Rows: thickness ratio (larger = thicker); Columns: angle');

    hasData = false;
    for hm = 1:numel(heatmapMetrics)
        metric = heatmapMetrics{hm};
        mat = data.(field).(metric);
        if ~all(isnan(mat), 'all')
            hasData = true;
        end
        ax = nexttile(tl);
        imagesc(thetasDeg, trVals, mat);
        set(ax, 'YDir', 'normal');
        xlabel(ax, 'Lamellar angle (deg)');
        ylabel(ax, 'Thickness ratio');
        title(ax, sprintf('%s %s', fieldLabel, metric), 'Interpreter', 'tex');
        colorbar(ax);
    end
    if hasData
        exportgraphics(fig, fullfile(outDir, sprintf('heatmap_%s.png', lower(field))), 'Resolution', 300);
    end
    close(fig);
end

    function makeLines(field, vsRatio, fileName)
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, 1, numel(metricNames), 'TileSpacing', 'compact', 'Padding', 'compact');
        xLabelText = ternary(vsRatio, 'Thickness ratio', 'Lamellar angle (deg)');
        fieldLabel = format_field_label(field);
        title(tl, sprintf('%s aggregated (%s)', fieldLabel, ternary(vsRatio, 'vary ratio', 'vary angle')), ...
            'FontWeight', 'bold', 'Interpreter', 'tex');
        subtitle(tl, sprintf('Metrics across r \\in [%s], \\theta \\in [%s]', ...
            strjoin(ratioLabels, ', '), strjoin(angleLabels, ', ')));

        for m = 1:numel(metricNames)
            metric = metricNames{m};
            mat = data.(field).(metric);
            ax = nexttile(tl);
            hold(ax, 'on');
            if vsRatio
                for jTh = 1:nTheta
                    plot(ax, trVals, mat(:, jTh), '-o', 'DisplayName', sprintf('\\theta = %s', angleLabels{jTh}));
                end
            else
                for iTr = 1:nTr
                    plot(ax, thetasDeg, mat(iTr, :), '-o', 'DisplayName', ratioLabels{iTr});
                end
            end
            grid(ax, 'on');
            xlabel(ax, xLabelText);
            ylabel(ax, metric);
            title(ax, sprintf('%s %s', fieldLabel, metric), 'Interpreter', 'tex');
            legend(ax, 'Location', 'bestoutside');
        end

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        close(fig);
    end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function label = format_ratio_value(val)
if abs(round(val) - val) < 1e-9
    label = sprintf('%d', round(val));
else
    label = sprintf('%.2g', val);
    label = strrep(label, ' ', '');
end
end
