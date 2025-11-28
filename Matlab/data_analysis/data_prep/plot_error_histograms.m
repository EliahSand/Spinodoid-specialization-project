function plot_error_histograms(errorTable, metrics, trLabel, thetaDeg, outDir)
%PLOT_ERROR_HISTOGRAMS Plot absolute and relative error histograms.

arguments
    errorTable table
    metrics struct
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

ensure_dir(outDir);

% Strip down to the most telling field for curvature: U3 only.
fieldsList = {'U3'};
modelTag = sprintf('%s | \\theta = %d°', trLabel, thetaDeg);

maxMAEField = fieldsList{1};
maxRelField = fieldsList{1};

makeFigure(false, 'hist_abs_errors.png', maxMAEField);
makeFigure(true, 'hist_rel_errors.png', maxRelField);

    function makeFigure(isRelative, fileName, highlightField)
        name = fieldsList{1};
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        if isRelative
            title(tl, sprintf('Relative error histogram — %s', modelTag), 'FontWeight', 'bold');
        else
            title(tl, sprintf('Absolute error histogram — %s', modelTag), 'FontWeight', 'bold');
        end
        subtitle(tl, sprintf('Highlighted: %s', highlightField));

        ax = nexttile(tl);

        if isRelative
            data = abs(errorTable.(sprintf('rel_%s', name)));
        else
            data = abs(errorTable.(sprintf('d_%s', name)));
        end
        data(data == 0) = eps;

        histogram(ax, data, 'NumBins', 40, 'Normalization', 'probability');
        if ~isRelative
            set(ax, 'XScale', 'log');
        end
        grid(ax, 'on');

        mark = "";
        if strcmp(name, highlightField)
            mark = " ★";
        end
        m = metrics.(name);
        title(ax, sprintf('%s%s (MAE=%.3g, MaxRel=%.3g)', name, mark, m.MAE, m.MaxRel));
        if isRelative
            xlabel(ax, '|(Shell - Solid)/(Solid)|');
        else
            xlabel(ax, '|Shell - Solid|');
        end
        ylabel(ax, 'Probability');

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
