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
baseParams.bandwidth = 0.125;
baseParams.nModes = 4000;
baseParams.solid_frac = 0.50;
baseParams.coneDeg = [30 0 0];
baseParams.rngSeed = 3;
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

% ----------------------------------------------------------------
% Local-defect sweep (cracks + holes, no GRF).
% Builds cfg.local_cases: 35 cases = 5 baselines + 15 cracks + 15 holes.
% ----------------------------------------------------------------
crackSev.length_m  = 25e-3;  % 62.5% of L (was 15e-3)
crackSev.theta_deg = 0;
crackSev.width_m   = 2e-3;   % wider so the band is clearly visible

holeSev.radius_m = 5e-3;     % radius 12.5% of L (was 3e-3)
holeSev.count    = 1;

% Canonical positions stored as fractions of L.
canonPos = struct( ...
    'name', {'center',    'quarter'}, ...
    'frac', {[0.50 0.50], [0.25 0.25]});

crackPosNames = {'center', 'quarter'};
holePosNames  = {'center', 'quarter'};

nLocalPerAngle = 1 + numel(crackPosNames) + numel(holePosNames);  % 7
nLocalTotal    = nAng * nLocalPerAngle;  % 35

local_cases = struct( ...
    'case_name',          cell(nLocalTotal, 1), ...
    'case_dir',           cell(nLocalTotal, 1), ...
    'angle_tag',          cell(nLocalTotal, 1), ...
    'lamellar_angle_deg', cell(nLocalTotal, 1), ...
    'defect_type',        cell(nLocalTotal, 1), ...
    'position_name',      cell(nLocalTotal, 1), ...
    'position_xy',        cell(nLocalTotal, 1), ...
    'severity',           cell(nLocalTotal, 1), ...
    'local_defects',      cell(nLocalTotal, 1));

k = 0;
for ai = 1:nAng
    ang    = lamellarAngles(ai);
    angTag = sprintf('lam%03d', ang);

    % Baseline (no defect)
    k = k + 1;
    local_cases(k).case_name          = [angTag '_baseline'];
    local_cases(k).case_dir           = fullfile(resultsRoot, angTag, 'baseline');
    local_cases(k).angle_tag          = angTag;
    local_cases(k).lamellar_angle_deg = ang;
    local_cases(k).defect_type        = 'baseline';
    local_cases(k).position_name      = 'none';
    local_cases(k).position_xy        = [];
    local_cases(k).severity           = struct();
    local_cases(k).local_defects      = struct([]);

    % Cracks at 3 positions
    for pi = 1:numel(crackPosNames)
        pName = crackPosNames{pi};
        pIdx  = find(strcmp({canonPos.name}, pName), 1);
        pXY   = canonPos(pIdx).frac * baseParams.L;

        ds          = struct();
        ds.type     = 'crack';
        ds.position = pXY;
        ds.length_m  = crackSev.length_m;
        ds.theta_deg = crackSev.theta_deg;
        ds.width_m   = crackSev.width_m;

        k = k + 1;
        local_cases(k).case_name          = [angTag '_crack_' pName];
        local_cases(k).case_dir           = fullfile(resultsRoot, angTag, ['crack_' pName]);
        local_cases(k).angle_tag          = angTag;
        local_cases(k).lamellar_angle_deg = ang;
        local_cases(k).defect_type        = 'crack';
        local_cases(k).position_name      = pName;
        local_cases(k).position_xy        = pXY;
        local_cases(k).severity           = crackSev;
        local_cases(k).local_defects      = ds;
    end

    % Holes at 3 positions
    for pi = 1:numel(holePosNames)
        pName = holePosNames{pi};
        pIdx  = find(strcmp({canonPos.name}, pName), 1);
        pXY   = canonPos(pIdx).frac * baseParams.L;

        ds          = struct();
        ds.type     = 'hole';
        ds.position = pXY;
        ds.radius_m  = holeSev.radius_m;
        ds.count     = holeSev.count;

        k = k + 1;
        local_cases(k).case_name          = [angTag '_hole_' pName];
        local_cases(k).case_dir           = fullfile(resultsRoot, angTag, ['hole_' pName]);
        local_cases(k).angle_tag          = angTag;
        local_cases(k).lamellar_angle_deg = ang;
        local_cases(k).defect_type        = 'hole';
        local_cases(k).position_name      = pName;
        local_cases(k).position_xy        = pXY;
        local_cases(k).severity           = holeSev;
        local_cases(k).local_defects      = ds;
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
cfg.local_cases     = local_cases;
cfg.baseParams      = baseParams;
end
