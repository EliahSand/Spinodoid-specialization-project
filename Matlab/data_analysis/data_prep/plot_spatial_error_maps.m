function plot_spatial_error_maps(solid, shell, errorTable, trLabel, thetaDeg, outDir)
%PLOT_SPATIAL_ERROR_MAPS Plot spatial maps for selected fields (solid/shell/error).

arguments
    solid table
    shell table
    errorTable table
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
end

ensure_dir(outDir);

coordsX = solid.X;
coordsY = solid.Z; % use X-Z plane consistently for all maps
modelTag = sprintf('%s | \\theta = %d°', trLabel, thetaDeg);

displacementFields = {'U1', 'U2', 'U3'};
stressFields = {'S11', 'S22', 'SMises'};

makeMapFigure(displacementFields, 'Displacement fields', 'map_displacements.png');
makeMapFigure(stressFields, 'Stress fields', 'map_stresses.png');

    function makeMapFigure(fieldList, figLabel, fileName)
        nRows = numel(fieldList);
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, nRows, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('%s — %s', figLabel, modelTag), 'FontWeight', 'bold');
        subtitle(tl, 'Columns: solid | shell | |error|');

        for i = 1:nRows
            fname = fieldList{i};
            solidVals = solid.(fname);
            shellVals = shell.(fname);
            errVals = abs(errorTable.(sprintf('d_%s', fname)));

            climVal = max(abs([solidVals; shellVals]), [], 'omitnan');
            if climVal == 0
                climVal = 1;
            end
            errMax = max(errVals, [], 'omitnan');
            if errMax == 0
                errMax = 1;
            end

            plotPanel(solidVals, [-climVal, climVal], sprintf('%s (solid)', fname), tl, 'parula');
            plotPanel(shellVals, [-climVal, climVal], sprintf('%s (shell)', fname), tl, 'parula');
            plotPanel(errVals, [0, errMax], sprintf('|d%s|', fname), tl, 'hot');
        end

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        close(fig);
    end

    function plotPanel(values, cLimits, ttl, layoutHandle, cmap)
        ax = nexttile(layoutHandle);
        scatter(ax, coordsX, coordsY, 18, values, 'filled', ...
            'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.85);
        axis(ax, 'equal');
        axis(ax, 'tight');
        grid(ax, 'on');
        colormap(ax, cmap);
        caxis(ax, cLimits);
        colorbar(ax);
        title(ax, ttl);
        xlabel(ax, 'X');
        ylabel(ax, 'Z');
    end
end
