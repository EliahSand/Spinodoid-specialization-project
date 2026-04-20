function cfg = defect_prediction_config()
%DEFECT_PREDICTION_CONFIG Shared settings for the defect prediction workflow.

scriptDir = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(scriptDir);
repoRoot = fileparts(matlabRoot);
resultsRoot = fullfile(scriptDir, 'results');

defectFractions = 0.0:0.1:1.0;
lamellarAngles  = [0 30 45 60 90];

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

% Build flat case list
nAng  = numel(lamellarAngles);
nFrac = numel(defectFractions);
cases = struct( ...
    'case_name',         cell(nAng * nFrac, 1), ...
    'case_dir',          cell(nAng * nFrac, 1), ...
    'angle_tag',         cell(nAng * nFrac, 1), ...
    'defect_tag',        cell(nAng * nFrac, 1), ...
    'lamellar_angle_deg', cell(nAng * nFrac, 1), ...
    'defect_fraction',   cell(nAng * nFrac, 1));

k = 0;
for ai = 1:nAng
    ang = lamellarAngles(ai);
    angTag = sprintf('lam%03d', ang);
    for fi = 1:nFrac
        frac = defectFractions(fi);
        defTag = sprintf('def%03d', round(100 * frac));
        k = k + 1;
        cases(k).angle_tag          = angTag;
        cases(k).defect_tag         = defTag;
        cases(k).case_name          = [angTag '_' defTag];
        cases(k).case_dir           = fullfile(resultsRoot, angTag, defTag);
        cases(k).lamellar_angle_deg = ang;
        cases(k).defect_fraction    = frac;
    end
end

cfg = struct();
cfg.root            = scriptDir;
cfg.resultsRoot     = resultsRoot;
cfg.matlabRoot      = matlabRoot;
cfg.repoRoot        = repoRoot;
cfg.stagingRoot     = fullfile(resultsRoot, '.staging');
cfg.defectFractions = defectFractions;
cfg.lamellarAngles  = lamellarAngles;
cfg.cases           = cases;
cfg.baseParams      = baseParams;
end
