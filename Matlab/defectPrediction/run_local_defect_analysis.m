% RUN_LOCAL_DEFECT_ANALYSIS  Rank lamellar angles by resistance to local defects.
%
%   Reads FEA results for the 35-case local-defect sweep and produces:
%     local_defect_summary.csv  — per-case metrics
%     angle_ranking.csv         — per-angle resistance scores and rank
%     angle_ranking_bar.png     — grouped bar chart (angles x defect types)
%     position_effect.png       — max|U3| per angle per position per type
%     FINDINGS_local.md         — auto-written text report
%
% How to run:
%   run('Matlab/defectPrediction/run_local_defect_analysis.m')

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

cfg = defect_prediction_config();
addpath(fullfile(cfg.matlabRoot, 'data_analysis', 'data_prep'));

outDir = fullfile(cfg.resultsRoot, 'analysis', 'local_defects');
ensure_dir(outDir);

cases  = cfg.local_cases;
angles = cfg.lamellarAngles;
nCases = numel(cases);

% -----------------------------------------------------------------
% 1. Load per-case metrics
% -----------------------------------------------------------------
template = struct( ...
    'case_name',          '', ...
    'lamellar_angle_deg', 0,  ...
    'defect_type',        '', ...
    'position_name',      '', ...
    'max_abs_u3_mm',      NaN, ...
    'p99_smises_pa',      NaN);
rows = repmat(template, 0, 1);

for k = 1:nCases
    c       = cases(k);
    csvPath = fullfile(c.case_dir, 'FEA', 'midplane_results.csv');
    if ~isfile(csvPath)
        fprintf('[SKIP] %s — no FEA results yet\n', c.case_name);
        continue;
    end

    tbl = load_spinodal_model(c.case_dir, 'solid');
    U3  = double(tbl.U3);
    Sm  = double(tbl.S_Mises);

    r                    = template;
    r.case_name          = c.case_name;
    r.lamellar_angle_deg = c.lamellar_angle_deg;
    r.defect_type        = c.defect_type;
    r.position_name      = c.position_name;
    r.max_abs_u3_mm      = max(abs(U3)) * 1e3;
    r.p99_smises_pa      = prctile(Sm, 99);
    rows(end+1) = r; %#ok<SAGROW>
    fprintf('[OK] %s\n', c.case_name);
end

if isempty(rows)
    fprintf('No FEA results found. Run FEA first.\n');
    return;
end
fprintf('Loaded %d / %d cases.\n\n', numel(rows), nCases);

% Write summary CSV.
summaryPath = fullfile(outDir, 'local_defect_summary.csv');
T = struct2table(rmfield(rows, {}));
writetable(T, summaryPath);
fprintf('Wrote %s\n', summaryPath);

% -----------------------------------------------------------------
% 2. Per-angle resistance scores
%    R(angle, type) = mean over positions of (max|U3|_defect / max|U3|_baseline)
% -----------------------------------------------------------------
defectTypes = {'crack', 'hole'};
nAng  = numel(angles);
nType = numel(defectTypes);

rankTemplate = struct( ...
    'lamellar_angle_deg', 0, ...
    'R_crack',  NaN, ...
    'R_hole',   NaN, ...
    'R_mean',   NaN, ...
    'rank',     NaN);
angle_metrics = repmat(rankTemplate, nAng, 1);

for ai = 1:nAng
    ang     = angles(ai);
    angMask = abs([rows.lamellar_angle_deg] - ang) < 0.5;
    angRows = rows(angMask);

    baseIdx = strcmp({angRows.defect_type}, 'baseline');
    if ~any(baseIdx)
        fprintf('[WARN] No baseline for angle %d — skipping.\n', ang);
        angle_metrics(ai).lamellar_angle_deg = ang;
        continue;
    end
    u3_base = angRows(baseIdx).max_abs_u3_mm;

    angle_metrics(ai).lamellar_angle_deg = ang;

    for ti = 1:nType
        dtype   = defectTypes{ti};
        defMask = strcmp({angRows.defect_type}, dtype);
        defRows = angRows(defMask);
        if isempty(defRows), continue; end
        R_vals = [defRows.max_abs_u3_mm] / max(u3_base, 1e-12);
        R_mean = mean(R_vals);
        angle_metrics(ai).(['R_' dtype]) = R_mean;
    end

    R_crack = angle_metrics(ai).R_crack;
    R_hole  = angle_metrics(ai).R_hole;
    valid   = [R_crack, R_hole];
    valid   = valid(~isnan(valid));
    angle_metrics(ai).R_mean = mean(valid);
end

% Rank by R_mean (ascending = most resistant first).
R_vals = [angle_metrics.R_mean];
[~, ord] = sort(R_vals, 'ascend', 'MissingPlacement', 'last');
for ri = 1:numel(ord)
    angle_metrics(ord(ri)).rank = ri;
end

rankPath = fullfile(outDir, 'angle_ranking.csv');
Tr = struct2table(angle_metrics);
writetable(Tr, rankPath);
fprintf('Wrote %s\n', rankPath);

% -----------------------------------------------------------------
% 3. Plots
% -----------------------------------------------------------------
angLabels = arrayfun(@(a) sprintf('%d°', a), angles, 'UniformOutput', false);
colors    = lines(nAng);

% — Grouped bar: R_crack and R_hole per angle —
R_crack_vals = [angle_metrics.R_crack];
R_hole_vals  = [angle_metrics.R_hole];
barData = [R_crack_vals(:), R_hole_vals(:)];

fig1 = figure('Visible','off');
b = bar(barData);
b(1).DisplayName = 'crack';
b(2).DisplayName = 'hole';
set(gca, 'XTickLabel', angLabels);
xlabel('Lamellar angle');
ylabel('Normalised max|U3| (lower = more resistant)');
title('Angle resistance to local defects');
legend('Location','northwest');
yline(1, '--k', 'Baseline', 'LabelVerticalAlignment','bottom');
saveas(fig1, fullfile(outDir, 'angle_ranking_bar.png'));
close(fig1);
fprintf('Wrote angle_ranking_bar.png\n');

% — Position effect: max|U3| per position per angle, one subplot per type —
fig2 = figure('Visible','off','Position',[0 0 1000 420]);
for ti = 1:nType
    dtype = defectTypes{ti};
    subplot(1, nType, ti);
    hold on;
    for ai = 1:nAng
        ang     = angles(ai);
        angMask = abs([rows.lamellar_angle_deg] - ang) < 0.5;
        defMask = strcmp({rows.defect_type}, dtype);
        sel     = rows(angMask & defMask);
        if isempty(sel), continue; end
        posNames = {sel.position_name};
        u3vals   = [sel.max_abs_u3_mm];
        [posNames, ord] = sort(posNames);
        u3vals = u3vals(ord);
        plot(1:numel(posNames), u3vals, '-o', 'Color', colors(ai,:), ...
            'DisplayName', angLabels{ai});
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', get_position_names(dtype, cfg));
    xlabel('Position');
    ylabel('max|U3| (mm)');
    title(['Position effect — ' dtype]);
    if ti == 1, legend('Location','best'); end
    hold off;
end
saveas(fig2, fullfile(outDir, 'position_effect.png'));
close(fig2);
fprintf('Wrote position_effect.png\n');

% -----------------------------------------------------------------
% 4. FINDINGS_local.md
% -----------------------------------------------------------------
write_findings(outDir, angle_metrics, angles, defectTypes, rows);
fprintf('Wrote FINDINGS_local.md\n\n');

% -----------------------------------------------------------------------
function names = get_position_names(dtype, cfg)
cases = cfg.local_cases;
typeMask = strcmp({cases.defect_type}, dtype) & ...
           abs([cases.lamellar_angle_deg] - cfg.lamellarAngles(1)) < 0.5;
sel = cases(typeMask);
[names, ord] = sort({sel.position_name});
end

function write_findings(outDir, angle_metrics, angles, defectTypes, rows)
fpath = fullfile(outDir, 'FINDINGS_local.md');
fid   = fopen(fpath, 'w');
if fid < 0, warning('Cannot write FINDINGS_local.md'); return; end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Local-Defect Resistance Findings\n\n');
fprintf(fid, 'Lamellar angles tested: %s degrees.\n\n', ...
    strjoin(arrayfun(@num2str, angles, 'UniformOutput', false), ', '));

for ti = 1:numel(defectTypes)
    dtype = defectTypes{ti};
    Rfield = ['R_' dtype];
    Rvals  = [angle_metrics.(Rfield)];
    valid  = ~isnan(Rvals);
    if ~any(valid), continue; end

    [~, bestIdx]  = min(Rvals(valid));
    [~, worstIdx] = max(Rvals(valid));
    validAngles   = angles(valid);
    bestAng       = validAngles(bestIdx);
    worstAng      = validAngles(worstIdx);
    bestR         = min(Rvals(valid));
    worstR        = max(Rvals(valid));

    fprintf(fid, '## %s\n\n', dtype);
    fprintf(fid, ...
        'Most resistant angle: **%d°** (mean normalised U3 = %.3f relative to baseline).\n', ...
        bestAng, bestR);
    fprintf(fid, ...
        'Least resistant angle: **%d°** (mean normalised U3 = %.3f).\n\n', ...
        worstAng, worstR);
end

% Overall ranking
fprintf(fid, '## Overall ranking (ascending = most resistant)\n\n');
[~, ord] = sort([angle_metrics.rank]);
for ri = 1:numel(ord)
    am = angle_metrics(ord(ri));
    fprintf(fid, '%d. %d° — R_crack=%.3f, R_hole=%.3f, R_mean=%.3f\n', ...
        am.rank, am.lamellar_angle_deg, am.R_crack, am.R_hole, am.R_mean);
end

% Note on position dominance
allR = [];
for ti = 1:numel(defectTypes)
    dtype = defectTypes{ti};
    for ai = 1:numel(angles)
        ang   = angles(ai);
        aMask = abs([rows.lamellar_angle_deg] - ang) < 0.5; %#ok<NODEF>
        dMask = strcmp({rows.defect_type}, dtype);
        sel   = rows(aMask & dMask);
        if ~isempty(sel)
            allR = [allR, std([sel.max_abs_u3_mm]) / max(mean([sel.max_abs_u3_mm]), 1e-12)]; %#ok<AGROW>
        end
    end
end
if ~isempty(allR)
    meanPosCoV = mean(allR);
    fprintf(fid, '\nMean within-angle CoV across positions: %.3f ', meanPosCoV);
    if meanPosCoV > 0.1
        fprintf(fid, '(position has a notable effect on response).\n');
    else
        fprintf(fid, '(position has minor effect; angle dominates).\n');
    end
end
end
