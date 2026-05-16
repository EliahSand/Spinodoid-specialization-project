% preview_v3_raster_boundary
%
% Build the V3 raster-augmented graph for ONE sample in memory (no save, no
% dataset folder touched) and plot the boundary classification + edges.
% Lets you sanity-check changes to collect_raster_boundary_candidates before
% regenerating the full dataset.
%
% Edit INP_PATH (or leave '' to auto-pick the first sheet_shell.inp under
% data/raw/samples) and DETAIL_LEVEL, then run.

scriptDir = fileparts(mfilename('fullpath'));
gnnRoot   = fileparts(scriptDir);
repoRoot  = fileparts(fileparts(gnnRoot));
addpath(fullfile(gnnRoot, 'helpers'));
addpath(fullfile(gnnRoot, 'datasetCreation'));
addpath(fullfile(gnnRoot, 'datasetCreation', 'gnn_graph'));
gnnPrepRoot = fullfile(repoRoot, 'Matlab', 'gnn_prep_spinodal');
addpath(gnnPrepRoot);
addpath(genpath(fullfile(gnnPrepRoot, 'src')));

INP_PATH     = '';      % '' = auto-pick first available
DETAIL_LEVEL = 0.3;

if isempty(INP_PATH)
    rawRoot = fullfile(gnnRoot, 'data', 'raw', 'samples');
    found = dir(fullfile(rawRoot, '**', 'sheet_shell.inp'));
    if isempty(found)
        error('No sheet_shell.inp found under %s', rawRoot);
    end
    INP_PATH = fullfile(found(1).folder, found(1).name);
end
fprintf('Previewing: %s\n', INP_PATH);

inpData = read_abaqus_inp(INP_PATH);
fullGraph = build_full_reference_graph_gnn(inpData, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'AutoDetectElset', true, ...
    'PatternPriority', {'spinodal', 'top'});
structuralGraph = extract_structural_graph(fullGraph, ...
    'DetailLevel', DETAIL_LEVEL, ...
    'MinIslandNodes', 1);
dense_data = local_rasterize(fullGraph, 128);

xy  = structuralGraph.node_coords;
xLo = min(xy(:,1)); xHi = max(xy(:,1));
yLo = min(xy(:,2)); yHi = max(xy(:,2));

[newXY, newCls, newParent, relIdx, relCls] = ...
    local_compute_v3_boundary_augmentation(xy, dense_data);

classIdxCenter = ones(size(xy, 1), 1);   % interior default
for ri = 1:numel(relIdx)
    cIdx = relIdx(ri); newC = relCls(ri);
    if classIdxCenter(cIdx) == 1 || newC < classIdxCenter(cIdx)
        classIdxCenter(cIdx) = newC;
    end
end

Nc = size(xy, 1);
Nb = size(newXY, 1);
xy_aug   = [xy; newXY];
classAll = [classIdxCenter; newCls];

ei_old = double(structuralGraph.edge_index);
newEdges = zeros(2, Nb);
for b = 1:Nb
    newEdges(:, b) = [newParent(b); Nc + b];
end
ei_aug = [ei_old, newEdges];

fprintf('Centerline nodes: %d   Added boundary nodes: %d (total %d)\n', Nc, Nb, Nc + Nb);
fprintf('Per-class totals: interior=%d clamped=%d loaded=%d free=%d\n', ...
    nnz(classAll==1), nnz(classAll==2), nnz(classAll==3), nnz(classAll==4));

classNames = {'interior', 'clamped\_left', 'loaded\_right', 'free\_top\_bot'};
classCols  = [0.75 0.75 0.75;
              0.85 0.10 0.10;
              0.10 0.45 0.85;
              0.20 0.65 0.20];

figure('Color', 'w', 'Position', [100 100 760 700]); hold on;

% raster occupancy as faint background
occ = logical(dense_data.raster(:, :, 1));
[Hg, Wg] = size(occ);
xGrid = linspace(xLo, xHi, Wg);
yGrid = linspace(yLo, yHi, Hg);
imagesc(xGrid, yGrid, double(occ), 'AlphaData', 0.15);
colormap(gca, gray);

% edges (excluded from legend)
edgeX = nan(3, size(ei_aug, 2));
edgeY = nan(3, size(ei_aug, 2));
edgeX(1, :) = xy_aug(ei_aug(1, :), 1).';
edgeX(2, :) = xy_aug(ei_aug(2, :), 1).';
edgeY(1, :) = xy_aug(ei_aug(1, :), 2).';
edgeY(2, :) = xy_aug(ei_aug(2, :), 2).';
plot(edgeX(:), edgeY(:), '-', 'Color', [0 0 0 0.25], 'LineWidth', 0.4, ...
    'HandleVisibility', 'off');

% nodes coloured by class
hScatter = gobjects(1, 4);
for c = 1:4
    sel = classAll == c;
    hScatter(c) = scatter(xy_aug(sel, 1), xy_aug(sel, 2), 22, classCols(c, :), 'filled', ...
        'DisplayName', sprintf('%s (n=%d)', classNames{c}, nnz(sel)));
end

axis equal tight; box on;
xlabel('x'); ylabel('y');
[~, runName] = fileparts(fileparts(INP_PATH));
title(sprintf('V3 preview (detail=%.2f) — %s', DETAIL_LEVEL, strrep(runName, '_', '\_')));
legend(hScatter, 'Location', 'eastoutside');
set(gca, 'YDir', 'normal');

function dense_data = local_rasterize(fullGraph, gridSize)
xy = double(fullGraph.node_coords(:, 1:2));
boundary = logical(fullGraph.boundary_mask(:));
xLo = min(xy(:, 1)); xHi = max(xy(:, 1));
yLo = min(xy(:, 2)); yHi = max(xy(:, 2));
cols = 1 + round((xy(:, 1) - xLo) / max(xHi - xLo, eps) * (gridSize - 1));
rows = 1 + round((xy(:, 2) - yLo) / max(yHi - yLo, eps) * (gridSize - 1));
cols = min(gridSize, max(1, cols));
rows = min(gridSize, max(1, rows));
occupancy = false(gridSize, gridSize);
boundaryMask = false(gridSize, gridSize);
idx = sub2ind([gridSize, gridSize], rows, cols);
occupancy(idx) = true;
boundaryMask(idx(boundary)) = true;
if isfield(fullGraph, 'edges_local') && ~isempty(fullGraph.edges_local)
    edges = double(fullGraph.edges_local);
    for e = 1:size(edges, 1)
        s = edges(e, 1); t = edges(e, 2);
        nSteps = max(abs(rows(s) - rows(t)), abs(cols(s) - cols(t))) + 1;
        rr = round(linspace(rows(s), rows(t), nSteps));
        cc = round(linspace(cols(s), cols(t), nSteps));
        edgeIdx = sub2ind([gridSize, gridSize], rr, cc);
        occupancy(edgeIdx) = true;
    end
end
dense_data = struct();
dense_data.raster = uint8(cat(3, occupancy, boundaryMask));
dense_data.xy_bounds = [xLo, xHi, yLo, yHi];
dense_data.grid_size = gridSize;
end

function [newXY, newCls, newParent, relabelIdx, relabelCls] = ...
        local_compute_v3_boundary_augmentation(xy, dense_data)
% Mirror of compute_v3_boundary_augmentation in step3_batch_structural_graph_from_inp.m.
occ = logical(dense_data.raster(:, :, 1));
[H, W] = size(occ);
xLo = dense_data.xy_bounds(1); xHi = dense_data.xy_bounds(2);
yLo = dense_data.xy_bounds(3); yHi = dense_data.xy_bounds(4);
dxPix = (xHi - xLo) / max(W - 1, 1);
dyPix = (yHi - yLo) / max(H - 1, 1);
relabelMaxDist = 1.5 * max(dxPix, dyPix);

walls = struct('xy', {}, 'cls', {});
rL = find(occ(:, 1));
if ~isempty(rL)
    walls(end+1).xy = [repmat(xLo, numel(rL), 1), yLo + (rL - 1) * dyPix];
    walls(end).cls  = 2;
end
rR = find(occ(:, W));
if ~isempty(rR)
    walls(end+1).xy = [repmat(xHi, numel(rR), 1), yLo + (rR - 1) * dyPix];
    walls(end).cls  = 3;
end
cB = find(occ(1, :));
if ~isempty(cB)
    walls(end+1).xy = [xLo + (cB(:) - 1) * dxPix, repmat(yLo, numel(cB), 1)];
    walls(end).cls  = 4;
end
cT = find(occ(H, :));
if ~isempty(cT)
    walls(end+1).xy = [xLo + (cT(:) - 1) * dxPix, repmat(yHi, numel(cT), 1)];
    walls(end).cls  = 4;
end

newXY      = zeros(0, 2);
newCls     = zeros(0, 1);
newParent  = zeros(0, 1);
relabelIdx = zeros(0, 1);
relabelCls = zeros(0, 1);

for w = 1:numel(walls)
    pxy = walls(w).xy;
    cls = walls(w).cls;
    k   = size(pxy, 1);
    parent = zeros(k, 1);
    for p = 1:k
        d2 = sum((xy - pxy(p, :)).^2, 2);
        [~, parent(p)] = min(d2);
    end
    [uParent, ~, grp] = unique(parent);
    for g = 1:numel(uParent)
        sel = find(grp == g);
        cIdx = uParent(g);
        if numel(sel) >= 2
            newXY(end+1, :)     = mean(pxy(sel, :), 1); %#ok<AGROW>
            newCls(end+1, 1)    = cls;                  %#ok<AGROW>
            newParent(end+1, 1) = cIdx;                 %#ok<AGROW>
        else
            pPos = pxy(sel, :);
            if norm(xy(cIdx, :) - pPos) <= relabelMaxDist
                relabelIdx(end+1, 1) = cIdx;            %#ok<AGROW>
                relabelCls(end+1, 1) = cls;             %#ok<AGROW>
            else
                newXY(end+1, :)     = pPos;             %#ok<AGROW>
                newCls(end+1, 1)    = cls;              %#ok<AGROW>
                newParent(end+1, 1) = cIdx;             %#ok<AGROW>
            end
        end
    end
end
end
