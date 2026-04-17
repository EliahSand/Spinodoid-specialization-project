% RUN_COMPARE_DEFECT_DEFORMATION Compare defect-case solid deformations visually.
%
% How to run:
%   run('Matlab/defectPrediction/run_compare_defect_deformation.m')
%
% Run this after:
%   abaqus python Matlab/defectPrediction/run_batch_defect_abaqus.py
%
% Outputs are written under:
%   Matlab/defectPrediction/results/analysis/deformation

cfg = defect_prediction_config();
dataPrepDir = fullfile(cfg.matlabRoot, 'data_analysis', 'data_prep');
addpath(dataPrepDir);

outDir = fullfile(cfg.resultsRoot, 'analysis', 'deformation');
ensure_dir(outDir);

profiles = struct( ...
    'case_name', {}, ...
    'defect_fraction', {}, ...
    'Y0_mm', {}, ...
    'Z0_mm', {}, ...
    'Ydef_mm', {}, ...
    'Zdef_mm', {}, ...
    'U3_mm', {}, ...
    'max_abs_u3_mm', {}, ...
    'rms_u3_mm', {}, ...
    'peak_to_peak_u3_mm', {});

for i = 1:numel(cfg.caseDirs)
    caseDir = cfg.caseDirs{i};
    caseName = cfg.caseNames{i};
    defectFrac = cfg.defectFractions(i);
    csvPath = fullfile(caseDir, 'FEA', 'midplane_results.csv');
    if ~isfile(csvPath)
        fprintf('[SKIP] %s: missing %s\n', caseName, csvPath);
        continue;
    end

    tbl = load_spinodal_model(caseDir, 'solid');
    profiles(end+1) = build_profile(tbl, caseName, defectFrac); %#ok<SAGROW>
end

if isempty(profiles)
    error('run_compare_defect_deformation:NoResults', ...
        'No solid midplane_results.csv files were found under %s.', cfg.resultsRoot);
end

profiles = sort_profiles_by_defect(profiles);
baseline = get_baseline_profile(profiles);

summary = struct2table(rmfield(profiles, {'Y0_mm', 'Z0_mm', 'Ydef_mm', 'Zdef_mm', 'U3_mm'}));
summary = add_percent_deviation_columns(summary, baseline);
writetable(summary, fullfile(outDir, 'defect_deformation_summary.csv'));
save(fullfile(outDir, 'defect_deformation_summary.mat'), 'summary', 'profiles', 'baseline');

plot_u3_overlay(profiles, outDir);
plot_absolute_deviation_overlay(profiles, baseline, outDir);

fprintf('Saved defect deformation comparison outputs to: %s\n', outDir);

function profiles = sort_profiles_by_defect(profiles)
[~, order] = sort([profiles.defect_fraction]);
profiles = profiles(order);
end

function baseline = get_baseline_profile(profiles)
idx = find(abs([profiles.defect_fraction]) < 1e-12, 1, 'first');
if isempty(idx)
    error('run_compare_defect_deformation:MissingBaseline', ...
        'No def000 baseline case was found in the loaded profiles.');
end
baseline = profiles(idx);
end

function summary = add_percent_deviation_columns(summary, baseline)
baselineMax = baseline.max_abs_u3_mm;
baselineRms = baseline.rms_u3_mm;
baselineP2P = baseline.peak_to_peak_u3_mm;

summary.max_abs_u3_pct_vs_def000 = percent_delta(summary.max_abs_u3_mm, baselineMax);
summary.rms_u3_pct_vs_def000 = percent_delta(summary.rms_u3_mm, baselineRms);
summary.peak_to_peak_u3_pct_vs_def000 = percent_delta(summary.peak_to_peak_u3_mm, baselineP2P);
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
U3 = double(tbl.U3(:));

[Y0u, Z0u] = average_by_reference(Y0, Z0);
[~, Ydefu] = average_by_reference(Y0, Ydef);
[~, Zdefu] = average_by_reference(Y0, Zdef);
[~, U3u] = average_by_reference(Y0, U3);

[Y0u, order] = sort(Y0u);
Z0u = Z0u(order);
Ydefu = Ydefu(order);
Zdefu = Zdefu(order);
U3u = U3u(order);

profile = struct();
profile.case_name = caseName;
profile.defect_fraction = defectFrac;
profile.Y0_mm = 1e3 * Y0u;
profile.Z0_mm = 1e3 * Z0u;
profile.Ydef_mm = 1e3 * Ydefu;
profile.Zdef_mm = 1e3 * Zdefu;
profile.U3_mm = 1e3 * U3u;
profile.max_abs_u3_mm = max(abs(profile.U3_mm));
profile.rms_u3_mm = sqrt(mean(profile.U3_mm.^2));
profile.peak_to_peak_u3_mm = max(profile.U3_mm) - min(profile.U3_mm);
end

function [xUnique, yMean] = average_by_reference(x, y)
[xUnique, ~, ic] = unique(x(:), 'sorted');
yMean = accumarray(ic, y(:), [], @mean);
end

function plot_u3_overlay(profiles, outDir)
nProfiles = numel(profiles);
colors = turbo(nProfiles);

fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
for i = 1:nProfiles
    plot(ax, profiles(i).Y0_mm, profiles(i).U3_mm, ...
        'LineWidth', 1.8, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s (%.1f)', ...
        profiles(i).case_name, profiles(i).defect_fraction));
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]');
ylabel(ax, 'U_3 [mm]');
title(ax, 'Out-of-plane deformation vs defect fraction', 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside');
exportgraphics(fig, fullfile(outDir, 'u3_overlay.png'), 'Resolution', 300);
close(fig);
end

function plot_absolute_deviation_overlay(profiles, baseline, outDir)
nProfiles = numel(profiles);
colors = turbo(nProfiles);
yQuery = baseline.Y0_mm(:);
baselineU3 = baseline.U3_mm(:);

fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
for i = 1:nProfiles
    u3Interp = interp1(profiles(i).Y0_mm, profiles(i).U3_mm, ...
        yQuery, 'linear', 'extrap');
    absDeviation = u3Interp - baselineU3;
    plot(ax, yQuery, absDeviation, ...
        'LineWidth', 1.8, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('%s (%.1f)', ...
        profiles(i).case_name, profiles(i).defect_fraction));
end
grid(ax, 'on');
xlabel(ax, 'Y_0 [mm]');
ylabel(ax, '\DeltaU_3 relative to def000 [mm]');
title(ax, 'Absolute deformation deviation from def000', 'FontWeight', 'bold');
legend(ax, 'Location', 'eastoutside');
exportgraphics(fig, fullfile(outDir, 'u3_absolute_deviation_from_def000.png'), 'Resolution', 300);
close(fig);
end
