% run_single_curvature_analysis
% Minimal single-run analysis entry point.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

% Example run folder (relative to repository root):
runFolder = fullfile('Matlab', 'results', 'sheets', ...
    'lamellar', 'sheetCone_tr50_ang045_lamellar_N128_1x1');

if ~isfolder(runFolder)
    error(['Edit runFolder in run_single_curvature_analysis.m first. ' ...
        'Current path not found: %s'], runFolder);
end

compare_single_run(runFolder);
