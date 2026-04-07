function plot_scatter_fields(solid, shell, trLabel, thetaDeg, outDir)
%PLOT_SCATTER_FIELDS Create solid-vs-shell overlay scatter plots.

arguments
    solid table
    shell table
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

ensure_dir(outDir);

displacementFields = {'U1', 'U2', 'U3'};
normalStressFields = {'S_11', 'S_22', 'S_33', 'S_Mises'};
shearStressFields = {'S_12', 'S_13', 'S_23'};

ratioLabel = format_ratio_label(trLabel);
thetaLabel = format_theta_label(thetaDeg);
modelTag = sprintf('%s | %s', ratioLabel, thetaLabel);

% Displacements: one wide figure
makeScatterFigure(displacementFields, 'Overlay: displacements', ...
    'overlay_displacements.png', 3, 'compact', 'loose', 12);

% Stresses: split normal and shear for readability
makeScatterFigure(normalStressFields, 'Overlay: stresses (normal + Mises)', ...
    'overlay_stresses_normal.png', 2, 'loose', 'compact', 20);
makeScatterFigure(shearStressFields, 'Overlay: stresses (shear)', ...
    'overlay_stresses_shear.png', 2, 'loose', 'compact', 20);

    function makeScatterFigure(fieldList, figLabel, fileName, nCols, spacing, padding, markerSize)
        n = numel(fieldList);
        nCols = min(nCols, n);
        nRows = ceil(n / nCols);
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', spacing, 'Padding', padding);
        title(tl, sprintf('%s — %s', figLabel, modelTag), 'FontWeight', 'bold', ...
            'FontSize', 12, 'Interpreter', 'tex');
        subtitle(tl, 'Shell vs solid');

        for i = 1:n
            name = fieldList{i};
            if ~ismember(name, solid.Properties.VariableNames) || ~ismember(name, shell.Properties.VariableNames)
                continue;
            end

            ax = nexttile(tl);
            x = double(solid.(name));
            y = double(shell.(name));
            d = y - x;

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

            % Compact quality indicators for each overlay plot
            mae = mean(abs(d), 'omitnan');
            rVal = NaN;
            if std(x) > 0 && std(y) > 0
                cc = corrcoef(x, y);
                if size(cc, 1) >= 2
                    rVal = cc(1, 2);
                end
            end

            fieldLabel = format_field_label(name);
            title(ax, {fieldLabel; sprintf('R=%.3f  MAE=%.3g', rVal, mae)}, ...
                'FontSize', 10, 'Interpreter', 'tex');
            xlabel(ax, sprintf('%s (solid)', fieldLabel), 'Interpreter', 'tex');
            ylabel(ax, sprintf('%s (shell)', fieldLabel), 'Interpreter', 'tex');
        end

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        close(fig);
    end
end

function out = format_theta_label(thetaDeg)
if isnan(thetaDeg)
    out = '\theta = N/A';
else
    out = sprintf('\\theta = %d°', round(thetaDeg));
end
end
