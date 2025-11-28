function summary = run_all_analysis()
%RUN_ALL_ANALYSIS Compare solid vs shell models across all ratios/angles.
%   summary = run_all_analysis() drives loading, metrics, plotting, PCA, and
%   aggregated reporting for the lamellar sheet data set.

tic;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

trLabels = {'tr50', 'tr100', 'tr200'};
thetasDeg = [0, 30, 45, 60, 90];

doPCA = license('test', 'statistics_toolbox') && exist('pca', 'file') == 2;

resultsRoot = fullfile(scriptDir, '..', 'results');
analysisRoot = fullfile(resultsRoot, 'analysis');
summaryRoot = fullfile(resultsRoot, 'summary');
ensure_dir(analysisRoot);
ensure_dir(summaryRoot);

overlayData = struct('U1', struct(), 'U2', struct(), 'U3', struct());

summaryTemplate = struct('trLabel', '', 'thetaDeg', NaN, 'metrics', [], ...
    'status', 'pending', 'dispAccept', false, 'stressAccept', false, ...
    'dispMeanRel', NaN, 'stressMeanRel', NaN, 'overallScore', NaN);
summary(numel(trLabels), numel(thetasDeg)) = summaryTemplate; %#ok<AGROW>

dispThresh = struct('MaxRel', 0.15, 'R', 0.95);
stressThresh = struct('MaxRel', 0.25, 'R', 0.90);

fprintf('--- Running solid vs shell comparison across %d ratios x %d angles ---\n', ...
    numel(trLabels), numel(thetasDeg));

for iTr = 1:numel(trLabels)
    for jTh = 1:numel(thetasDeg)
        trLabel = trLabels{iTr};
        theta = thetasDeg(jTh);
        fprintf('\n[%s | \\theta = %d°] Loading data...\n', trLabel, theta);
        try
            [solid, shell] = load_solid_shell_pair(trLabel, theta);
            [metrics, errorTable] = compare_solid_shell(solid, shell);
            metricsTbl = metrics_to_table(metrics);

            outputDir = fullfile(analysisRoot, trLabel, sprintf('ang%03d', theta));
            ensure_dir(outputDir);
            writetable(metricsTbl, fullfile(outputDir, 'metrics_summary.csv'));
            writetable(errorTable, fullfile(outputDir, 'node_errors.csv'));
            save(fullfile(outputDir, 'metrics.mat'), 'metrics');
            save(fullfile(outputDir, 'error_table.mat'), 'errorTable');

            plot_scatter_fields(solid, shell, metrics, trLabel, theta, outputDir);
            plot_error_histograms(errorTable, metrics, trLabel, theta, outputDir);
            % Spatial maps proved uninformative for these midsections; rely on
            % curvature profiles and overlays instead.
            plot_curvature_profiles(solid, shell, metrics, trLabel, theta, outputDir);
            if doPCA
                analyze_error_patterns(solid, shell, trLabel, theta, outputDir);
            end

            dispFields = {'U1', 'U2', 'U3'};
            stressFields = {'S11', 'S22', 'SMises'};
            dispAccept = all(cellfun(@(f) metrics.(f).MaxRel < dispThresh.MaxRel ...
                && metrics.(f).R > dispThresh.R, dispFields));
            stressAccept = all(cellfun(@(f) metrics.(f).MaxRel < stressThresh.MaxRel ...
                && metrics.(f).R > stressThresh.R, stressFields));

            summary(iTr, jTh).trLabel = trLabel;
            summary(iTr, jTh).thetaDeg = theta;
            summary(iTr, jTh).metrics = metrics;
            summary(iTr, jTh).status = 'ok';
            summary(iTr, jTh).dispAccept = dispAccept;
            summary(iTr, jTh).stressAccept = stressAccept;
            summary(iTr, jTh).dispMeanRel = mean(cellfun(@(f) metrics.(f).MeanRel, dispFields));
            summary(iTr, jTh).stressMeanRel = mean(cellfun(@(f) metrics.(f).MeanRel, stressFields));
            summary(iTr, jTh).overallScore = mean([summary(iTr, jTh).dispMeanRel, summary(iTr, jTh).stressMeanRel]);

            % Store displacement data for overlay plots per thickness ratio
            fieldsToStore = {'U1','U2','U3'};
            for fIdx = 1:numel(fieldsToStore)
                fname = fieldsToStore{fIdx};
                if ~isfield(overlayData.(fname), trLabel)
                    overlayData.(fname).(trLabel) = struct('thetaDeg', [], 'solid', {{}}, 'shell', {{}}); %#ok<CCAT>
                end
                overlayData.(fname).(trLabel).thetaDeg(end+1) = theta;
                overlayData.(fname).(trLabel).solid{end+1} = double(solid.(fname));
                overlayData.(fname).(trLabel).shell{end+1} = double(shell.(fname));
            end

            fprintf('Component-level metrics:\n');
            print_component_summary(metrics, {'U1','U2','U3','S11','S22','SMises'});
            fprintf('  Displacements: %s (MaxRel < %.2f, R > %.2f)\n', ...
                ternary(dispAccept, 'ACCEPTABLE', 'NOT acceptable'), dispThresh.MaxRel, dispThresh.R);
            fprintf('  Stresses:      %s (MaxRel < %.2f, R > %.2f)\n', ...
                ternary(stressAccept, 'ACCEPTABLE', 'NOT acceptable'), stressThresh.MaxRel, stressThresh.R);
        catch ME
            summary(iTr, jTh).trLabel = trLabel;
            summary(iTr, jTh).thetaDeg = theta;
            summary(iTr, jTh).status = sprintf('failed: %s', ME.message);
            stackMsg = "";
            if ~isempty(ME.stack)
                s = ME.stack(1);
                [~, fname, ext] = fileparts(s.file);
                stackMsg = sprintf(' (%s%s:%d)', fname, ext, s.line);
            end
            warning('Failed for %s @ %d°: %s%s', trLabel, theta, ...
                ME.getReport('extended', 'hyperlinks', 'off'), stackMsg);
        end
    end
end

overviewTbl = build_overview_table(summary);
writetable(overviewTbl, fullfile(summaryRoot, 'summary_overview.csv'));
save(fullfile(summaryRoot, 'comparison_summary.mat'), 'summary', 'trLabels', 'thetasDeg');
plot_aggregated_metrics(summary, trLabels, thetasDeg, summaryRoot);
plot_overlay_series(overlayData, summaryRoot);

toc;

report_global_summary(overviewTbl);
end

function print_component_summary(metrics, fields)
for i = 1:numel(fields)
    name = fields{i};
    m = metrics.(name);
    fprintf('  %-6s MAE=%8.3g | MaxRel=%7.3g | R=%6.3f\n', ...
        name, m.MAE, m.MaxRel, m.R);
end
end

function overviewTbl = build_overview_table(summary)
rows = [];
for i = 1:numel(summary)
    s = summary(i);
    rows = [rows; {s.trLabel, s.thetaDeg, s.status, s.dispMeanRel, s.stressMeanRel, ...
        s.dispAccept, s.stressAccept, s.overallScore}]; %#ok<AGROW>
end
overviewTbl = cell2table(rows, 'VariableNames', ...
    {'trLabel','thetaDeg','status','dispMeanRel','stressMeanRel','dispAccept','stressAccept','overallScore'});
end

function report_global_summary(overviewTbl)
okIdx = strcmp(overviewTbl.status, 'ok');
if ~any(okIdx)
    fprintf('\nNo successful model comparisons to summarize.\n');
    return;
end

validTbl = overviewTbl(okIdx, :);
[~, bestDispIdx] = min(validTbl.dispMeanRel);
[~, bestStressIdx] = min(validTbl.stressMeanRel);

bestDisp = validTbl(bestDispIdx, :);
bestStress = validTbl(bestStressIdx, :);

fprintf('\n=== Global summary ===\n');
fprintf('Best displacement agreement: %s @ %d° (mean rel = %.3g)\n', ...
    bestDisp.trLabel{1}, bestDisp.thetaDeg, bestDisp.dispMeanRel);
fprintf('Best stress agreement:        %s @ %d° (mean rel = %.3g)\n', ...
    bestStress.trLabel{1}, bestStress.thetaDeg, bestStress.stressMeanRel);

highDisp = validTbl(validTbl.dispMeanRel > 0.5, :);
highStress = validTbl(validTbl.stressMeanRel > 0.5, :);
if ~isempty(highDisp)
    fprintf('Combinations with large displacement deviations (>0.5 mean rel):\n');
    disp(highDisp(:, {'trLabel','thetaDeg','dispMeanRel'}));
end
if ~isempty(highStress)
    fprintf('Combinations with large stress deviations (>0.5 mean rel):\n');
    disp(highStress(:, {'trLabel','thetaDeg','stressMeanRel'}));
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
