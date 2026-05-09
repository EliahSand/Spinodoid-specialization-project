% visualize_junction_spline_graph.m
%
% Demo visualization of a single junction-spline (v5) graph.
% Shows the JL graph overlaid on the 128x128 occupancy mask.
%
% Usage:
%   visualize_junction_spline_graph           % picks a random sample
%   sampleName = 'tr30_ang45_...'; visualize_junction_spline_graph

scriptDir  = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(scriptDir);
gnnRoot    = fullfile(matlabRoot, 'GNN');
samplesDir = fullfile(gnnRoot, 'data', 'dataset_junction_spline', 'samples');

if ~exist('sampleName', 'var') || isempty(sampleName)
    items = dir(samplesDir);
    items = items([items.isdir] & ~startsWith({items.name}, '.'));
    if isempty(items)
        error('No samples found in %s\nRun step3_batch_junction_spline first.', samplesDir);
    end
    idx = randi(numel(items));
    sampleName = items(idx).name;
    fprintf('Picked random sample: %s\n', sampleName);
end

d = load(fullfile(samplesDir, sampleName, 'sample.mat'));
jl   = d.jlGraph;
dens = d.dense_data;

occ = double(dens.raster(:,:,1));   % 128x128 occupancy
bnd = double(dens.raster(:,:,2));   % 128x128 boundary

xy   = jl.node_coords;             % nJL x 2  (physical)
ef   = jl.edge_features;           % E x 4  [length, mean_r, cos2t, sin2t]
ei   = jl.edges_local;             % E x 2
isBnd = jl.node_features(:,4);     % nJL x 1

% Map physical coords -> raster pixel coords for background overlay
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));
G   = 128;
col = 1 + (xy(:,1)-xLo)/max(xHi-xLo,eps)*(G-1);
row = 1 + (xy(:,2)-yLo)/max(yHi-yLo,eps)*(G-1);

% Edge colormap: phys_length
if ~isempty(ef)
    eLen = ef(:,1);
    eRad = ef(:,2);
    theta = atan2(ef(:,4), ef(:,3)) / 2;   % recover theta from (cos2t,sin2t)
else
    eLen = []; eRad = []; theta = [];
end

figure('Name', sprintf('Junction-spline graph: %s', sampleName), ...
    'Position', [50 50 1300 500]);

colorEdgesBy = {'Length', 'Mean radius', 'Orientation'};
edgeData     = {eLen, eRad, theta * 180/pi};
cmaps        = {parula, cool, hsv};

for panel = 1:3
    subplot(1,3,panel);
    % Background: fixed RGB image so the axes colormap (used by colorbar)
    % does not corrupt the background. Light = material, dark = void.
    bg = repmat(0.15 + 0.65 * occ, 1, 1, 3);   % H x W x 3, values in [0.15, 0.80]
    image(bg); axis image; hold on;
    set(gca,'YDir','normal');

    if ~isempty(ei)
        cData = edgeData{panel};
        cLo = min(cData); cHi = max(cData);
        if cHi - cLo < eps, cHi = cLo + 1; end
        cNorm = (cData - cLo) / (cHi - cLo);
        cmap  = cmaps{panel};
        cIdx  = max(1, min(256, round(cNorm * 255) + 1));
        cRGB  = cmap(cIdx, :);

        for e = 1:size(ei,1)
            u = ei(e,1); v = ei(e,2);
            if u == v
                plot(col(u), row(u), 'o', 'Color', cRGB(e,:), 'MarkerSize', 10, 'LineWidth', 2);
            else
                plot([col(u), col(v)], [row(u), row(v)], '-', ...
                    'Color', cRGB(e,:), 'LineWidth', 2);
            end
        end
        colormap(gca, cmap);
        clim([cLo, cHi]);
        cb = colorbar; cb.Label.String = colorEdgesBy{panel};
    end

    % Nodes: junctions = filled circle, endpoints = diamond
    % degree computed from edge list
    nodeDeg = accumarray(ei(:), 1, [jl.num_nodes, 1]);
    isEnd   = nodeDeg == 1 | nodeDeg == 0;

    scatter(col(~isEnd & ~isBnd), row(~isEnd & ~isBnd), 40, [0.2 0.6 1], 'filled', 'MarkerEdgeColor','k','LineWidth',0.5);
    scatter(col(isEnd  & ~isBnd), row(isEnd  & ~isBnd), 50, [1 0.5 0],   'd', 'filled', 'MarkerEdgeColor','k','LineWidth',0.5);
    scatter(col(logical(isBnd)),  row(logical(isBnd)),  50, [0.9 0.2 0.2],'s', 'filled', 'MarkerEdgeColor','k','LineWidth',0.5);

    title(sprintf('Edges by %s', colorEdgesBy{panel}), 'FontSize', 10);
    axis off;
end

sgtitle(sprintf('%s  |  %d nodes  |  %d edges', ...
    strrep(sampleName,'_','\_'), jl.num_nodes, jl.num_edges), 'FontSize', 11);

% --- Edge feature distributions ---
if ~isempty(ef) && size(ef,1) > 1
    figure('Name','Edge feature distributions','Position',[50 600 900 300]);
    fnames = {'phys\_length', 'mean\_radius', 'cos2\theta', 'sin2\theta'};
    for k = 1:4
        subplot(1,4,k);
        histogram(ef(:,k), 30, 'FaceColor', [0.3 0.6 0.9]);
        title(fnames{k}); xlabel('value'); ylabel('count');
        grid on; box off;
    end
    sgtitle(sprintf('Edge feature distributions  (%d edges)', size(ef,1)));
end

fprintf('Nodes: %d  |  Edges: %d\n', jl.num_nodes, jl.num_edges);
fprintf('Edge length:  min=%.3g  max=%.3g  mean=%.3g\n', min(eLen), max(eLen), mean(eLen));
fprintf('Mean radius:  min=%.3g  max=%.3g  mean=%.3g\n', min(eRad), max(eRad), mean(eRad));
