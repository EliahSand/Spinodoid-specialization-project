% RUN_COMPARE_LOCAL_DEFECTS  Thesis-quality comparison of local-defect deformation.
%
%   Produces, for the 35-case local-defect sweep:
%     analysis/local_defects_compare/
%       lam{ang}/u3_overlay.png              — 7 cases overlaid (baseline + 6 defects)
%       lam{ang}/u3_deviation.png            — defect minus baseline (6 curves)
%       lam{ang}/geometry_strip.png          — top-down spin-layer for all 7 cases
%       by_defect/u3_overlay_{defect}.png    — same defect across 5 angles
%       insensitivity_heatmap.png            — max|U3|_defect / max|U3|_base (angle x case)
%       deviation_heatmap.png                — max|ΔU3| / max|U3|_base (angle x case)
%       angle_insensitivity_bar.png          — per-angle mean deviation score
%       local_defect_compare_summary.csv     — per-case scalar metrics
%       FINDINGS_local_compare.md            — narrative argument for thesis
%
% How to run (from repo root or from Matlab/):
%   run('Matlab/defectPrediction/run_compare_local_defects.m')
%
% Requires FEA results: abaqus python run_batch_defect_abaqus.py (local cases).

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

cfg = defect_prediction_config();
addpath(fullfile(cfg.matlabRoot, 'data_analysis', 'data_prep'));

outDir = fullfile(cfg.resultsRoot, 'analysis', 'local_defects_compare');
ensure_dir(outDir);

cases  = cfg.local_cases;
angles = cfg.lamellarAngles;

% -----------------------------------------------------------------
% 1. Load per-case profiles
% -----------------------------------------------------------------
profiles = struct( ...
    'case_name',          {}, ...
    'lamellar_angle_deg', {}, ...
    'angle_tag',          {}, ...
    'defect_type',        {}, ...
    'position_name',      {}, ...
    'defect_tag',         {}, ...
    'case_dir',           {}, ...
    'Y_mm',               {}, ...
    'U3_mm',              {}, ...
    'max_abs_u3_mm',      {}, ...
    'rms_u3_mm',          {}, ...
    'peak_to_peak_u3_mm', {});

for k = 1:numel(cases)
    c = cases(k);
    csvPath = fullfile(c.case_dir, 'FEA', 'midplane_results.csv');
    if ~isfile(csvPath)
        fprintf('[SKIP] %s — no FEA results\n', c.case_name);
        continue;
    end
    tbl = load_spinodal_model(c.case_dir, 'solid');
    [Yu, U3u] = average_by_reference(double(tbl.Y), double(tbl.U3));
    [Yu, ord] = sort(Yu);
    U3u = U3u(ord);

    p = struct();
    p.case_name          = c.case_name;
    p.lamellar_angle_deg = c.lamellar_angle_deg;
    p.angle_tag          = c.angle_tag;
    p.defect_type        = c.defect_type;
    p.position_name      = c.position_name;
    p.defect_tag         = compose_defect_tag(c.defect_type, c.position_name);
    p.case_dir           = c.case_dir;
    p.Y_mm               = 1e3 * Yu;
    p.U3_mm              = 1e3 * U3u;
    p.max_abs_u3_mm      = max(abs(p.U3_mm));
    p.rms_u3_mm          = sqrt(mean(p.U3_mm.^2));
    p.peak_to_peak_u3_mm = max(p.U3_mm) - min(p.U3_mm);
    profiles(end+1) = p; %#ok<SAGROW>
end

if isempty(profiles)
    error('run_compare_local_defects:NoResults', ...
        'No FEA results found under %s. Run the FEA batch first.', cfg.resultsRoot);
end
fprintf('Loaded %d / %d cases.\n', numel(profiles), numel(cases));

defectTagsAll = {'baseline', 'crack_center', 'crack_quarter', 'crack_mid_edge_x', ...
                 'hole_center',  'hole_quarter',  'hole_mid_edge_y'};
defectColors  = thesis_defect_colors();

% -----------------------------------------------------------------
% 2. Per-angle overlays + deviation + geometry strip
% -----------------------------------------------------------------
for ai = 1:numel(angles)
    ang     = angles(ai);
    angTag  = sprintf('lam%03d', ang);
    angDir  = fullfile(outDir, angTag);
    ensure_dir(angDir);

    angMask = abs([profiles.lamellar_angle_deg] - ang) < 0.5;
    angProfs = profiles(angMask);
    if isempty(angProfs)
        fprintf('[SKIP] %s — no profiles\n', angTag);
        continue;
    end

    baseIdx = find(strcmp({angProfs.defect_type}, 'baseline'), 1);
    if isempty(baseIdx)
        fprintf('[WARN] %s — no baseline; deviation skipped\n', angTag);
        baseline = [];
    else
        baseline = angProfs(baseIdx);
    end

    plot_u3_overlay_named(angProfs, defectTagsAll, defectColors, baseline, ...
        sprintf('%s — U_3 along midplane (x = L/2)', angTag), ...
        fullfile(angDir, 'u3_overlay.png'));

    if ~isempty(baseline)
        plot_deviation_named(angProfs, baseline, defectTagsAll, defectColors, ...
            sprintf('%s — \\DeltaU_3 vs baseline', angTag), ...
            fullfile(angDir, 'u3_deviation.png'));
    end

    plot_geometry_strip(angProfs, defectTagsAll, ...
        sprintf('%s — spin-layer top-down (per case)', angTag), ...
        fullfile(angDir, 'geometry_strip.png'));
end

% -----------------------------------------------------------------
% 3. Per-defect: same defect across all angles
% -----------------------------------------------------------------
byDefectDir = fullfile(outDir, 'by_defect');
ensure_dir(byDefectDir);
angleColors = lines(numel(angles));

for di = 1:numel(defectTagsAll)
    dtag = defectTagsAll{di};
    sel  = profiles(strcmp({profiles.defect_tag}, dtag));
    if numel(sel) < 2, continue; end

    [~, ord] = sort([sel.lamellar_angle_deg]);
    sel = sel(ord);

    fig = figure('Visible', 'off', 'Position', [0 0 900 420]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');
    for i = 1:numel(sel)
        cidx = find(angles == sel(i).lamellar_angle_deg, 1);
        plot(ax, sel(i).Y_mm, sel(i).U3_mm, ...
            'LineWidth', 1.8, 'Color', angleColors(cidx, :), ...
            'DisplayName', sprintf('%d°', sel(i).lamellar_angle_deg));
    end
    grid(ax, 'on');
    xlabel(ax, 'Y_0 [mm]'); ylabel(ax, 'U_3 [mm]');
    title(ax, sprintf('%s — U_3 across lamellar angles', dtag), ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    legend(ax, 'Location', 'eastoutside');
    exportgraphics(fig, fullfile(byDefectDir, sprintf('u3_overlay_%s.png', dtag)), 'Resolution', 200);
    close(fig);
end

% -----------------------------------------------------------------
% 4. Build (angle x defect) matrices and emit heatmaps + bar
% -----------------------------------------------------------------
defectCols = setdiff(defectTagsAll, {'baseline'}, 'stable');
nAng = numel(angles); nDef = numel(defectCols);

R_ins  = nan(nAng, nDef);   % max|U3|_defect / max|U3|_baseline   (1 = insensitive)
R_dev  = nan(nAng, nDef);   % max|ΔU3|        / max|U3|_baseline   (0 = insensitive)
RMS_R  = nan(nAng, nDef);   % rms ΔU3         / max|U3|_baseline

for ai = 1:nAng
    ang     = angles(ai);
    angMask = abs([profiles.lamellar_angle_deg] - ang) < 0.5;
    angProfs = profiles(angMask);
    bIdx = find(strcmp({angProfs.defect_type}, 'baseline'), 1);
    if isempty(bIdx), continue; end
    base = angProfs(bIdx);
    u3b_max = max(abs(base.U3_mm));

    for di = 1:nDef
        dtag = defectCols{di};
        m = strcmp({angProfs.defect_tag}, dtag);
        if ~any(m), continue; end
        p = angProfs(m);
        u3i = interp1(p.Y_mm, p.U3_mm, base.Y_mm, 'linear', 'extrap');
        dU3 = u3i - base.U3_mm;
        R_ins(ai, di) = max(abs(p.U3_mm)) / max(u3b_max, eps);
        R_dev(ai, di) = max(abs(dU3))     / max(u3b_max, eps);
        RMS_R(ai, di) = sqrt(mean(dU3.^2)) / max(u3b_max, eps);
    end
end

plot_heatmap(R_ins, angles, defectCols, ...
    'max |U_3|_{defect} / max |U_3|_{baseline}   (1 = no change)', ...
    fullfile(outDir, 'insensitivity_heatmap.png'), [0.7 1.3]);

devTop = max(R_dev(:), [], 'omitnan');
if isempty(devTop) || ~isfinite(devTop) || devTop <= 0
    devClim = [];
else
    devClim = [0 devTop * 1.05];
end
plot_heatmap(R_dev, angles, defectCols, ...
    'max |\DeltaU_3| / max |U_3|_{baseline}   (0 = no change)', ...
    fullfile(outDir, 'deviation_heatmap.png'), devClim);

% Per-angle bar: mean deviation across defects (lower = more insensitive).
meanDev = mean(R_dev, 2, 'omitnan');
fig = figure('Visible', 'off', 'Position', [0 0 720 420]);
b = bar(meanDev); %#ok<NASGU>
set(gca, 'XTickLabel', arrayfun(@(a) sprintf('%d°', a), angles, 'UniformOutput', false));
xlabel('Lamellar angle');
ylabel('mean_{defects}  max|\DeltaU_3| / max|U_3|_{baseline}');
title('Defect insensitivity by lamellar angle (lower = more insensitive)', 'FontWeight', 'bold');
grid on;
yline(0, '--k');
for i = 1:numel(meanDev)
    if ~isnan(meanDev(i))
        text(i, meanDev(i), sprintf(' %.2f', meanDev(i)), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
exportgraphics(fig, fullfile(outDir, 'angle_insensitivity_bar.png'), 'Resolution', 200);
close(fig);

% -----------------------------------------------------------------
% 5. Summary CSV + findings markdown
% -----------------------------------------------------------------
write_summary_csv(profiles, fullfile(outDir, 'local_defect_compare_summary.csv'));
write_findings(outDir, angles, defectCols, R_ins, R_dev, RMS_R, meanDev);

fprintf('Saved local-defect comparison outputs to: %s\n', outDir);

% =========================================================================
%  HELPERS
% =========================================================================

function tag = compose_defect_tag(defect_type, position_name)
if strcmpi(defect_type, 'baseline')
    tag = 'baseline';
else
    tag = sprintf('%s_%s', defect_type, position_name);
end
end

function colors = thesis_defect_colors()
% Visually distinct colors for the 7 case categories.
colors = struct( ...
    'baseline',         [0 0 0], ...
    'crack_center',     [0.84 0.10 0.10], ...
    'crack_quarter',    [0.93 0.43 0.17], ...
    'crack_mid_edge_x', [0.99 0.65 0.20], ...
    'hole_center',      [0.13 0.36 0.73], ...
    'hole_quarter',     [0.30 0.60 0.86], ...
    'hole_mid_edge_y',  [0.45 0.78 0.95]);
end

function [xUnique, yMean] = average_by_reference(x, y)
[xUnique, ~, ic] = unique(x(:), 'sorted');
yMean = accumarray(ic, y(:), [], @mean);
end

function plot_u3_overlay_named(profs, tagOrder, colorStruct, baseline, titleStr, outPath)
fig = figure('Visible', 'off', 'Position', [0 0 900 480]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
for ti = 1:numel(tagOrder)
    tag = tagOrder{ti};
    m = strcmp({profs.defect_tag}, tag);
    if ~any(m), continue; end
    p = profs(m);
    if strcmp(tag, 'baseline')
        plot(ax, p.Y_mm, p.U3_mm, 'Color', colorStruct.(tag), ...
            'LineWidth', 3.0, 'DisplayName', 'baseline');
    else
        plot(ax, p.Y_mm, p.U3_mm, 'Color', colorStruct.(tag), ...
            'LineWidth', 1.6, 'DisplayName', tag);
    end
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]'); ylabel(ax, 'U_3 [mm]');
title(ax, titleStr, 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
if ~isempty(baseline)
    yMax = max(abs(baseline.U3_mm));
    if yMax > 0
        text(ax, 0.02, 0.97, sprintf('max |U_3|_{base} = %.3f mm', yMax), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1]);
    end
end
exportgraphics(fig, outPath, 'Resolution', 200);
close(fig);
end

function plot_deviation_named(profs, baseline, tagOrder, colorStruct, titleStr, outPath)
fig = figure('Visible', 'off', 'Position', [0 0 900 480]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
yQ = baseline.Y_mm(:);
b  = baseline.U3_mm(:);
yMaxBase = max(abs(b));
for ti = 1:numel(tagOrder)
    tag = tagOrder{ti};
    if strcmp(tag, 'baseline'), continue; end
    m = strcmp({profs.defect_tag}, tag);
    if ~any(m), continue; end
    p = profs(m);
    u3i = interp1(p.Y_mm, p.U3_mm, yQ, 'linear', 'extrap');
    dU = u3i - b;
    plot(ax, yQ, dU, 'Color', colorStruct.(tag), ...
        'LineWidth', 1.8, 'DisplayName', sprintf('%s (max=%.2f%%)', ...
        tag, 100 * max(abs(dU)) / max(yMaxBase, eps)));
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]'); ylabel(ax, '\DeltaU_3 = U_3^{defect} - U_3^{base} [mm]');
title(ax, titleStr, 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
yline(ax, 0, '--k', 'HandleVisibility', 'off');
exportgraphics(fig, outPath, 'Resolution', 200);
close(fig);
end

function plot_geometry_strip(angProfs, tagOrder, titleStr, outPath)
% Render each case as a top-down image: light blue = spinodoid material,
% red = defect (material removed vs baseline), white = void.
baseIdx = find(strcmp({angProfs.defect_tag}, 'baseline'), 1);
baselineTD = [];
if ~isempty(baseIdx)
    matPath = fullfile(angProfs(baseIdx).case_dir, 'sheet.mat');
    if isfile(matPath)
        Sb = load(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness');
        baselineTD = compute_top_down(Sb);
    end
end

tiles = {};
labels = {};
outlinesPerTile = {};
for ti = 1:numel(tagOrder)
    tag = tagOrder{ti};
    m = strcmp({angProfs.defect_tag}, tag);
    if ~any(m), continue; end
    p = angProfs(m);
    matPath = fullfile(p.case_dir, 'sheet.mat');
    if ~isfile(matPath), continue; end
    S = load(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness');
    [td, Lx] = compute_top_down(S);
    if ~isempty(baselineTD) && isequal(size(baselineTD), size(td)) && ~strcmp(tag, 'baseline')
        defectMask = baselineTD & ~td;
    else
        defectMask = false(size(td));
    end
    rgb = compose_rgb(td, defectMask);

    outlines = {};
    metaPath = fullfile(p.case_dir, 'defect_case.json');
    if isfile(metaPath)
        meta = jsondecode(fileread(metaPath));
        if isfield(meta, 'local_defects') && ~isempty(meta.local_defects)
            outlines = defect_outlines(meta.local_defects);
        end
    end

    tiles{end+1} = struct('rgb', rgb, 'Lx', Lx); %#ok<AGROW>
    labels{end+1} = tag; %#ok<AGROW>
    outlinesPerTile{end+1} = outlines; %#ok<AGROW>
end
n = numel(tiles);
if n == 0, return; end

fig = figure('Visible', 'off', 'Position', [0 0 220 * n + 80, 280]);
tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:n
    ax = nexttile(tl);
    image(ax, [0 tiles{i}.Lx], [0 tiles{i}.Lx], tiles{i}.rgb);
    axis(ax, 'equal', 'tight');
    set(ax, 'YDir', 'normal');
    xlim(ax, [0 tiles{i}.Lx]); ylim(ax, [0 tiles{i}.Lx]);
    hold(ax, 'on');
    overlay_outlines(ax, outlinesPerTile{i}, tiles{i}.Lx);
    hold(ax, 'off');
    title(ax, labels{i}, 'Interpreter', 'none');
    if i == 1, ylabel(ax, 'y [mm]'); end
    xlabel(ax, 'x [mm]');
end
title(tl, titleStr, 'FontWeight', 'bold');
exportgraphics(fig, outPath, 'Resolution', 200);
close(fig);
end

function [topDown, Lx_mm] = compute_top_down(S)
% Returns topDown where row index ↔ y, column index ↔ x (matches the
% [row=y, col=x, z] convention used by apply_local_defects).
mask = S.sheetMask;
N    = size(mask, 1);
Svox = S.voxelSizeXY;
baseParams_L = 40e-3;
t_base = 2e-3;
tbV = max(3, round(t_base / (baseParams_L / N)));
spinLayer = mask(:,:,tbV+1:end);
topDown = squeeze(any(spinLayer, 3));
Lx_mm = N * Svox * 1e3;
end

function rgb = compose_rgb(solidMask, defectMask)
spinColor   = [0.60 0.78 0.95];   % light blue spinodoid
defectColor = [0.85 0.18 0.18];   % red defect
bgColor     = [1 1 1];

[H, W] = size(solidMask);
rgb = repmat(reshape(bgColor, 1, 1, 3), H, W);
solidIdx = solidMask & ~defectMask;
for c = 1:3
    chan = rgb(:,:,c);
    chan(solidIdx)  = spinColor(c);
    chan(defectMask) = defectColor(c);
    rgb(:,:,c) = chan;
end
end

function outlines = defect_outlines(localDefects)
% Build {[Nx2]} polylines (mm) for each defect's full intended footprint.
outlines = {};
specs = localDefects;
if iscell(specs), specs = [specs{:}]; end
for i = 1:numel(specs)
    d = specs(i);
    if ~isfield(d, 'type') || isempty(d.type), continue; end
    x0_mm = d.position(1) * 1e3;
    y0_mm = d.position(2) * 1e3;
    switch lower(d.type)
        case 'crack'
            L_mm = d.length_m * 1e3;
            W_mm = 0;
            if isfield(d, 'width_m'), W_mm = d.width_m * 1e3; end
            W_mm = max(W_mm, 0.2);
            th = deg2rad(d.theta_deg);
            cx = cos(th); sx = sin(th);
            hx = (L_mm/2) * [cx; sx];
            hw = (W_mm/2) * [-sx; cx];
            c = [hx + hw, -hx + hw, -hx - hw, hx - hw, hx + hw];
            outlines{end+1} = [x0_mm + c(1,:); y0_mm + c(2,:)]'; %#ok<AGROW>
        case 'hole'
            r_mm  = d.radius_m * 1e3;
            count = 1;
            if isfield(d, 'count'), count = d.count; end
            spiral = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1];
            th = linspace(0, 2*pi, 96);
            for ci = 1:count
                off = spiral(min(ci, size(spiral,1)), :) * 2.5 * r_mm;
                outlines{end+1} = [x0_mm + off(1) + r_mm * cos(th); ...
                                   y0_mm + off(2) + r_mm * sin(th)]'; %#ok<AGROW>
            end
    end
end
end

function overlay_outlines(ax, outlines, Lx)
if isempty(outlines), return; end
offsets = [0 0; Lx 0; -Lx 0; 0 Lx; 0 -Lx; Lx Lx; -Lx Lx; Lx -Lx; -Lx -Lx];
for i = 1:numel(outlines)
    P = outlines{i};
    for k = 1:size(offsets, 1)
        plot(ax, P(:,1) + offsets(k,1), P(:,2) + offsets(k,2), ...
            ':', 'Color', [0.85 0.18 0.18], 'LineWidth', 1.4);
    end
end
end

function plot_heatmap(M, angles, defectCols, cbLabel, outPath, climRange)
fig = figure('Visible', 'off', 'Position', [0 0 900 400]);
ax = axes(fig); %#ok<LAXES>
imagesc(ax, M);
if nargin >= 6 && ~isempty(climRange) && all(isfinite(climRange))
    caxis(ax, climRange);
end
colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = cbLabel; cb.Label.Interpreter = 'tex';
set(ax, 'XTick', 1:numel(defectCols), 'XTickLabel', defectCols, ...
        'YTick', 1:numel(angles),     'YTickLabel', ...
        arrayfun(@(a) sprintf('%d°', a), angles, 'UniformOutput', false));
xtickangle(ax, 30);
ylabel(ax, 'Lamellar angle');
title(ax, cbLabel, 'FontWeight', 'bold', 'Interpreter', 'tex');
% Annotate cells.
[na, nd] = size(M);
for i = 1:na
    for j = 1:nd
        if isnan(M(i,j)), continue; end
        text(ax, j, i, sprintf('%.2f', M(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 9, 'Color', 'w');
    end
end
exportgraphics(fig, outPath, 'Resolution', 200);
close(fig);
end

function write_summary_csv(profiles, outPath)
T = table( ...
    string({profiles.case_name}'), ...
    [profiles.lamellar_angle_deg]', ...
    string({profiles.defect_type}'), ...
    string({profiles.position_name}'), ...
    [profiles.max_abs_u3_mm]', ...
    [profiles.rms_u3_mm]', ...
    [profiles.peak_to_peak_u3_mm]', ...
    'VariableNames', {'case_name', 'lamellar_angle_deg', 'defect_type', ...
                      'position_name', 'max_abs_u3_mm', 'rms_u3_mm', 'peak_to_peak_u3_mm'});
writetable(T, outPath);
end

function write_findings(outDir, angles, defectCols, R_ins, R_dev, RMS_R, meanDev)
fpath = fullfile(outDir, 'FINDINGS_local_compare.md');
fid = fopen(fpath, 'w'); if fid < 0, return; end
c = onCleanup(@() fclose(fid));

fprintf(fid, '# Local-defect deformation comparison\n\n');
fprintf(fid, 'Per-case max |ΔU3| relative to that angle''s baseline max |U3|.\n');
fprintf(fid, 'Low values support the argument that the spinodoid shell is **defect-insensitive**.\n\n');

% Overall stats
allDev = R_dev(:); allDev = allDev(~isnan(allDev));
if ~isempty(allDev)
    fprintf(fid, '## Overall insensitivity\n\n');
    fprintf(fid, '- Mean  max|ΔU3| / max|U3|_base across all (angle, defect) pairs: **%.1f%%**\n', 100*mean(allDev));
    fprintf(fid, '- Median: **%.1f%%**\n', 100*median(allDev));
    fprintf(fid, '- 95th percentile: **%.1f%%**\n', 100*prctile(allDev, 95));
    fprintf(fid, '- Worst case: **%.1f%%**\n\n', 100*max(allDev));
end

% Best angle
[~, bestIdx]  = min(meanDev);
[~, worstIdx] = max(meanDev);
fprintf(fid, '## Ranking — most insensitive lamellar angle\n\n');
fprintf(fid, '| Rank | Angle | mean max|ΔU3| / max|U3|_base |\n|---|---|---|\n');
[~, ord] = sort(meanDev, 'ascend');
for ri = 1:numel(ord)
    fprintf(fid, '| %d | %d° | %.3f |\n', ri, angles(ord(ri)), meanDev(ord(ri)));
end
fprintf(fid, '\nMost insensitive: **%d°** (mean %.1f%%). Least insensitive: **%d°** (mean %.1f%%).\n\n', ...
    angles(bestIdx), 100*meanDev(bestIdx), angles(worstIdx), 100*meanDev(worstIdx));

% Per-defect worst angle
fprintf(fid, '## Worst angle per defect\n\n');
fprintf(fid, '| Defect | best angle | worst angle | best Δ%% | worst Δ%% |\n|---|---|---|---|---|\n');
for di = 1:numel(defectCols)
    col = R_dev(:, di);
    if all(isnan(col)), continue; end
    [bv, bi] = min(col, [], 'omitnan');
    [wv, wi] = max(col, [], 'omitnan');
    fprintf(fid, '| %s | %d° | %d° | %.1f%% | %.1f%% |\n', ...
        defectCols{di}, angles(bi), angles(wi), 100*bv, 100*wv);
end

fprintf(fid, '\n## Notes\n\n');
fprintf(fid, '- `insensitivity_heatmap.png` plots `max|U3|_def / max|U3|_base`; values close to 1.0 indicate the defect barely shifts peak deflection.\n');
fprintf(fid, '- `deviation_heatmap.png` plots `max|ΔU3| / max|U3|_base`; values close to 0.0 indicate the deflection profile is nearly preserved.\n');
fprintf(fid, '- `angle_insensitivity_bar.png` summarises the per-angle mean deviation; lower bars = stronger insensitivity argument.\n');
fprintf(fid, '- Per-angle figures in `lam{ang}/` show the actual U3 profile and ΔU3 — useful for showing a reader that defect curves visually overlap the baseline.\n');
fprintf(fid, '- `by_defect/u3_overlay_{tag}.png` lets you argue which angle handles a specific defect type best.\n');
end
