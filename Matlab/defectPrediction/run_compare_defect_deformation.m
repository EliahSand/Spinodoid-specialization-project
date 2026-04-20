% RUN_COMPARE_DEFECT_DEFORMATION Compare defect-case solid deformations visually.
%
% How to run:
%   run('Matlab/defectPrediction/run_compare_defect_deformation.m')
%
% Run this after:
%   abaqus python Matlab/defectPrediction/run_batch_defect_abaqus.py
%
% Outputs under Matlab/defectPrediction/results/analysis/deformation/:
%   defect_deformation_summary.csv / .mat       — all 55 cases
%   lam{ang}/u3_overlay.png                     — per-angle overlay (colored by fraction)
%   lam{ang}/u3_absolute_deviation_from_def000  — per-angle deviation from that angle's baseline
%   by_defect/u3_overlay_def{frac}.png          — cross-angle overlay per defect fraction
%   max_abs_u3_heatmap.png                      — angle × fraction heatmap

cfg = defect_prediction_config();
dataPrepDir = fullfile(cfg.matlabRoot, 'data_analysis', 'data_prep');
addpath(dataPrepDir);

outDir = fullfile(cfg.resultsRoot, 'analysis', 'deformation');
ensure_dir(outDir);

% --- Load all available profiles -----------------------------------------
profiles = struct( ...
    'case_name',          {}, ...
    'lamellar_angle_deg', {}, ...
    'defect_fraction',    {}, ...
    'angle_tag',          {}, ...
    'defect_tag',         {}, ...
    'Y0_mm', {}, 'Z0_mm', {}, 'Ydef_mm', {}, 'Zdef_mm', {}, 'U3_mm', {}, ...
    'max_abs_u3_mm', {}, 'rms_u3_mm', {}, 'peak_to_peak_u3_mm', {});

for i = 1:numel(cfg.cases)
    c = cfg.cases(i);
    csvPath = fullfile(c.case_dir, 'FEA', 'midplane_results.csv');
    if ~isfile(csvPath)
        fprintf('[SKIP] %s: missing %s\n', c.case_name, csvPath);
        continue;
    end

    tbl = load_spinodal_model(c.case_dir, 'solid');
    p = build_profile(tbl, c.case_name, c.defect_fraction);
    p.lamellar_angle_deg = c.lamellar_angle_deg;
    p.angle_tag          = c.angle_tag;
    p.defect_tag         = c.defect_tag;
    profiles(end+1) = p; %#ok<SAGROW>
end

if isempty(profiles)
    error('run_compare_defect_deformation:NoResults', ...
        'No solid midplane_results.csv files were found under %s.', cfg.resultsRoot);
end

% Sort by angle then fraction
[~, ord] = sortrows([[profiles.lamellar_angle_deg]' [profiles.defect_fraction]']);
profiles = profiles(ord);

% --- Summary table -------------------------------------------------------
scalar_fields = {'case_name', 'lamellar_angle_deg', 'defect_fraction', ...
    'max_abs_u3_mm', 'rms_u3_mm', 'peak_to_peak_u3_mm'};
summary = struct2table(rmfield(profiles, ...
    {'Y0_mm', 'Z0_mm', 'Ydef_mm', 'Zdef_mm', 'U3_mm', 'angle_tag', 'defect_tag'}));
summary = add_same_angle_pct_columns(summary, profiles);
writetable(summary, fullfile(outDir, 'defect_deformation_summary.csv'));
save(fullfile(outDir, 'defect_deformation_summary.mat'), 'summary', 'profiles');

% --- Per-angle plots ------------------------------------------------------
angles = cfg.lamellarAngles;
for ai = 1:numel(angles)
    ang = angles(ai);
    angTag = sprintf('lam%03d', ang);
    angProfiles = profiles(abs([profiles.lamellar_angle_deg] - ang) < 0.5);

    if isempty(angProfiles)
        fprintf('[SKIP] No profiles loaded for %s\n', angTag);
        continue;
    end

    angDir = fullfile(outDir, angTag);
    ensure_dir(angDir);

    % Baseline for this angle = def000
    baseIdx = find(abs([angProfiles.defect_fraction]) < 1e-12, 1, 'first');
    if isempty(baseIdx)
        fprintf('[WARN] No def000 baseline for %s; skipping deviation plot\n', angTag);
        baseline = [];
    else
        baseline = angProfiles(baseIdx);
    end

    fracs = [angProfiles.defect_fraction];
    colors = turbo(numel(angProfiles));

    plot_u3_overlay(angProfiles, fracs, colors, ...
        sprintf('%s — U_3 overlay', angTag), ...
        fullfile(angDir, 'u3_overlay.png'));

    if ~isempty(baseline)
        plot_deviation_overlay(angProfiles, baseline, fracs, colors, ...
            sprintf('%s — deviation from def000', angTag), ...
            fullfile(angDir, 'u3_absolute_deviation_from_def000.png'));
    end
end

% --- Cross-angle plots (one per defect fraction) -------------------------
byDefectDir = fullfile(outDir, 'by_defect');
ensure_dir(byDefectDir);
fracs = cfg.defectFractions;
angleColors = lines(numel(angles));

for fi = 1:numel(fracs)
    frac = fracs(fi);
    defTag = sprintf('def%03d', round(100 * frac));
    fracProfiles = profiles(abs([profiles.defect_fraction] - frac) < 1e-6);

    if isempty(fracProfiles)
        continue;
    end

    angVals = [fracProfiles.lamellar_angle_deg];
    % Color by angle index
    colorIdxs = arrayfun(@(a) find(angles == a, 1), angVals);
    pColors = angleColors(colorIdxs, :);

    plot_u3_overlay(fracProfiles, angVals, pColors, ...
        sprintf('def%03d — U_3 by lamellar angle', round(100 * frac)), ...
        fullfile(byDefectDir, sprintf('u3_overlay_%s.png', defTag)), ...
        'angle [deg]');
end

% --- Heatmap (angle x fraction) ------------------------------------------
plot_heatmap(profiles, angles, fracs, outDir);

fprintf('Saved defect deformation comparison outputs to: %s\n', outDir);

% =========================================================================
%  HELPERS
% =========================================================================

function summary = add_same_angle_pct_columns(summary, profiles)
pct_max  = nan(height(summary), 1);
pct_rms  = nan(height(summary), 1);
pct_p2p  = nan(height(summary), 1);

anglesPresent = unique([profiles.lamellar_angle_deg]);
for ai = 1:numel(anglesPresent)
    ang = anglesPresent(ai);
    mask = abs([profiles.lamellar_angle_deg] - ang) < 0.5;
    angProfiles = profiles(mask);
    baseIdx = find(abs([angProfiles.defect_fraction]) < 1e-12, 1);
    if isempty(baseIdx), continue; end
    base = angProfiles(baseIdx);

    rows = find(mask);
    pct_max(rows)  = percent_delta([profiles(mask).max_abs_u3_mm],  base.max_abs_u3_mm);
    pct_rms(rows)  = percent_delta([profiles(mask).rms_u3_mm],      base.rms_u3_mm);
    pct_p2p(rows)  = percent_delta([profiles(mask).peak_to_peak_u3_mm], base.peak_to_peak_u3_mm);
end

summary.max_abs_u3_pct_vs_def000_same_angle  = pct_max;
summary.rms_u3_pct_vs_def000_same_angle      = pct_rms;
summary.peak_to_peak_u3_pct_vs_def000_same_angle = pct_p2p;
end

function pct = percent_delta(values, baselineValue)
scale = max(abs(baselineValue), eps);
pct = 100 * (values - baselineValue) ./ scale;
end

function profile = build_profile(tbl, caseName, defectFrac)
Y0 = double(tbl.Y(:));
Z0 = double(tbl.Z(:));
Ydef = Y0 + double(tbl.U2(:));
Zdef = Z0 + double(tbl.U3(:));
U3   = double(tbl.U3(:));

[Y0u, Z0u]   = average_by_reference(Y0, Z0);
[~,   Ydefu] = average_by_reference(Y0, Ydef);
[~,   Zdefu] = average_by_reference(Y0, Zdef);
[~,   U3u]   = average_by_reference(Y0, U3);

[Y0u, order] = sort(Y0u);
Z0u   = Z0u(order);
Ydefu = Ydefu(order);
Zdefu = Zdefu(order);
U3u   = U3u(order);

profile = struct();
profile.case_name          = caseName;
profile.defect_fraction    = defectFrac;
profile.lamellar_angle_deg = 0;   % caller overwrites
profile.angle_tag          = '';
profile.defect_tag         = '';
profile.Y0_mm              = 1e3 * Y0u;
profile.Z0_mm              = 1e3 * Z0u;
profile.Ydef_mm            = 1e3 * Ydefu;
profile.Zdef_mm            = 1e3 * Zdefu;
profile.U3_mm              = 1e3 * U3u;
profile.max_abs_u3_mm      = max(abs(profile.U3_mm));
profile.rms_u3_mm          = sqrt(mean(profile.U3_mm.^2));
profile.peak_to_peak_u3_mm = max(profile.U3_mm) - min(profile.U3_mm);
end

function [xUnique, yMean] = average_by_reference(x, y)
[xUnique, ~, ic] = unique(x(:), 'sorted');
yMean = accumarray(ic, y(:), [], @mean);
end

function plot_u3_overlay(profs, colorValues, colors, titleStr, outPath, legendLabel)
if nargin < 6, legendLabel = 'defect frac'; end

fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
for i = 1:numel(profs)
    plot(ax, profs(i).Y0_mm, profs(i).U3_mm, ...
        'LineWidth', 1.8, 'Color', colors(i, :), ...
        'DisplayName', sprintf('%s=%.4g', legendLabel, colorValues(i)));
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]');
ylabel(ax, 'U_3 [mm]');
title(ax, titleStr, 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside');
exportgraphics(fig, outPath, 'Resolution', 300);
close(fig);
end

function plot_deviation_overlay(profs, baseline, fracVals, colors, titleStr, outPath)
fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
yQuery = baseline.Y0_mm(:);
baseU3 = baseline.U3_mm(:);
for i = 1:numel(profs)
    u3i = interp1(profs(i).Y0_mm, profs(i).U3_mm, yQuery, 'linear', 'extrap');
    plot(ax, yQuery, u3i - baseU3, ...
        'LineWidth', 1.8, 'Color', colors(i, :), ...
        'DisplayName', sprintf('frac=%.1f', fracVals(i)));
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]');
ylabel(ax, '\DeltaU_3 relative to def000 [mm]');
title(ax, titleStr, 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside');
exportgraphics(fig, outPath, 'Resolution', 300);
close(fig);
end

function plot_heatmap(profiles, angles, fracs, outDir)
nAng  = numel(angles);
nFrac = numel(fracs);
Z = nan(nAng, nFrac);

for ai = 1:nAng
    for fi = 1:nFrac
        mask = abs([profiles.lamellar_angle_deg] - angles(ai)) < 0.5 & ...
               abs([profiles.defect_fraction]    - fracs(fi))  < 1e-6;
        if any(mask)
            Z(ai, fi) = profiles(mask).max_abs_u3_mm;
        end
    end
end

fig = figure('Visible', 'off', 'Position', [0 0 900 400]);
ax = axes(fig); %#ok<LAXES>
imagesc(ax, fracs, angles, Z);
colormap(ax, 'turbo');
cb = colorbar(ax);
cb.Label.String = 'max |U_3| [mm]';
set(ax, 'YTick', angles, 'XTick', fracs);
xlabel(ax, 'defect fraction');
ylabel(ax, 'lamellar angle [deg]');
title(ax, 'max |U_3| — angle × defect fraction', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outDir, 'max_abs_u3_heatmap.png'), 'Resolution', 300);
close(fig);
end
