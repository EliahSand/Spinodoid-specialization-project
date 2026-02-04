function plot_scatter_fields(solid, shell, metrics, trLabel, thetaDeg, outDir)
%PLOT_SCATTER_FIELDS Create shell-vs-solid scatter plots for key components.

arguments
    solid table
    shell table
    metrics struct
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

ensure_dir(outDir);

displacementFields = {'U1', 'U2', 'U3'};
normalStressFields = {'S_11', 'S_22', 'S_33', 'S_Mises'};
shearStressFields = {'S_12', 'S_13', 'S_23'};

ratioLabel = format_ratio_label(trLabel);
modelTag = sprintf('%s | \\theta = %d°', ratioLabel, thetaDeg);

% Displacements: give extra padding so the title does not overlap axes
makeScatterFigure(displacementFields, 'Displacements', 'scatter_displacements.png', 3, 'compact', 'loose', 12);
% Stresses: split into two figures for readability
makeScatterFigure(normalStressFields, 'Stresses (normal & Mises)', 'scatter_stresses_normal.png', 2, 'loose', 'compact', 20);
makeScatterFigure(shearStressFields, 'Stresses (shear)', 'scatter_stresses_shear.png', 2, 'loose', 'compact', 20);

    function makeScatterFigure(fieldList, figLabel, fileName, nCols, spacing, padding, markerSize)
        n = numel(fieldList);
        nCols = min(nCols, n);
        nRows = ceil(n / nCols);
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', spacing, 'Padding', padding);
        title(tl, sprintf('%s — %s', figLabel, modelTag), 'FontWeight', 'bold', ...
            'FontSize', 12, 'Interpreter', 'tex');
        subtitle(tl, 'Shell vs. solid', 'FontSize', 10);

        for i = 1:n
            name = fieldList{i};
            ax = nexttile(tl);

            x = solid.(name);
            y = shell.(name);
            scatter(ax, x, y, markerSize, 'filled', 'MarkerFaceAlpha', 0.7);
            hold(ax, 'on');
            lims = [min([x; y], [], 'omitnan'), max([x; y], [], 'omitnan')];
            if diff(lims) == 0
                lims = lims + [-1, 1];
            else
                pad = 0.05 * diff(lims);
                lims = lims + [-pad, pad];
            end
            plot(ax, lims, lims, 'k--', 'LineWidth', 1.2);
            xlim(ax, lims);
            ylim(ax, lims);
            grid(ax, 'on');
            axis(ax, 'square');

            % Field title plus metrics on a second line
            m = metrics.(name);
            fieldLabel = format_field_label(name);
            title(ax, {fieldLabel; sprintf('R=%.3f  MAE=%.3g', m.R, m.MAE)}, ...
                'FontSize', 10, 'Interpreter', 'tex');
            xlabel(ax, sprintf('%s (solid)', fieldLabel), 'Interpreter', 'tex');
            ylabel(ax, sprintf('%s (shell)', fieldLabel), 'Interpreter', 'tex');

        end

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        close(fig);
    end
end
