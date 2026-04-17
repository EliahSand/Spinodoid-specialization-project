function cfg = defect_prediction_config()
%DEFECT_PREDICTION_CONFIG Shared settings for the defect prediction workflow.

scriptDir = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(scriptDir);
repoRoot = fileparts(matlabRoot);
resultsRoot = fullfile(scriptDir, 'results');

defectFractions = 0.0:0.1:1.0;
caseNames = arrayfun(@(f) sprintf('def%03d', round(100 * f)), ...
    defectFractions, 'UniformOutput', false);
caseDirs = cellfun(@(name) fullfile(resultsRoot, name), caseNames, ...
    'UniformOutput', false);

baseParams = struct();
baseParams.N = 128;
baseParams.L = 40e-3;
baseParams.lambda_vox = 25;
baseParams.bandwidth = 2;
baseParams.nModes = 4000;
baseParams.solid_frac = 0.50;
baseParams.coneDeg = [30 0 0];
baseParams.rngSeed = 1;
baseParams.t_spin = 1e-3;
baseParams.t_base = 2e-3;
baseParams.tilesXY = [1 1];
baseParams.slice_count = 8;
baseParams.lamellarAngleDeg = 0;
baseParams.remove_top_spin_frac = 0.0;

cfg = struct();
cfg.root = scriptDir;
cfg.resultsRoot = resultsRoot;
cfg.matlabRoot = matlabRoot;
cfg.repoRoot = repoRoot;
cfg.stagingRoot = fullfile(resultsRoot, '.staging');
cfg.defectFractions = defectFractions;
cfg.caseNames = caseNames;
cfg.caseDirs = caseDirs;
cfg.baseParams = baseParams;
end
