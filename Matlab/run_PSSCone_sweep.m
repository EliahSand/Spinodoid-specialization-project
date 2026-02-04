% RUN_PSSCONe_SWEEP  Batch runner for PSSCone with multiple lamellar angles
% and thickness ratios. Creates results under "results/sheets/lamellar/trXX".

function run_PSSCone_sweep()

tic

scriptDir = fileparts(mfilename('fullpath'));
resultsBase = fullfile(scriptDir, 'results', 'sheets', 'lamellar');
if ~exist(resultsBase, 'dir')
    mkdir(resultsBase);
end

angles = [0 30 45 60 90];           % lamellar angles in degrees
ratios = [0.5 1.0 2];               % t_spin / t_base ratios
baseThickness = 2e-3;               % meters

for r = ratios
    ratioLabel = sprintf('tr%02d', round(100 * r));
    ratioFolder = fullfile(resultsBase, ratioLabel);
    if ~exist(ratioFolder, 'dir')
        mkdir(ratioFolder);
    end

    params = struct();
    params.t_base = baseThickness;
    params.t_spin = r * baseThickness;
    params.resultsRoot = ratioFolder;
    params.tilesXY = [2 3];
    params.align_with_cube = true;

    % Parallel loop over angles (remember to start parallell pooling)
    parfor ia = 1:numel(angles)
        ang = angles(ia);
        p = params;                 %copy
        p.lamellarAngleDeg = ang;
        p.lamellarAngleDeg = ang;
        p.rngSeed = 1;
        fprintf('Running PSSCone: t_r=%.2f, angle=%.1f deg\n', r, ang);
        PSSCone(p);
    end

fprintf('---------------------------------------------- \n')
fprintf('FINISHED SWEEP\n')
toc
end
