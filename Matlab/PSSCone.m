% PSSCone - Spinodal two-layered sheet generator with controllable in-plane
% lamellar angle and tiling

function PSSCone(params)
tStart = tic;


if nargin < 1
    params = struct();
end

% ------------------------- Design parameters ---------------------------------------

cfg = struct();
cfg.N            = 2^7;          % grid size (NxNxN). Use powers of two for speed
cfg.L            = 40e-3;        % physical box length (mm)    (vanligvis 40e-3)

%   NB! Husk vox size is L/N

cfg.lambda_vox   = 25;          % target feature wavelength in voxels (~rib/ligament spacing)
cfg.bandwidth    = 2;           % relative shell thickness around target |k| (0.1–0.3)
cfg.nModes       = 4000;        % number of Fourier modes to sample (1k–10k typical)
cfg.solid_frac   = 0.50;        % volume fraction of SOLID after threshold (0..1)
cfg.coneDeg      = [30 0 0];    % cone half-angles about x,y,z [30 0 0] = lamellar 
cfg.rngSeed      = 1;           % reproducible
cfg.remove_top_spin_frac = 0.0;   % fraction of spinodal voxels to stochastically remove in the spinodoid layer

cfg.t_spin       = 1e-3;        %spinodal thickness         (vanligvis 1e-3)
cfg.t_base       = 2e-3;        %base thickness             (vanligvis 2e-3)
cfg.tilesXY      = [1 1];       %tiling for periodicity
cfg.slice_count  = 8;           %builds 2D spinodal patterm for shell by averaging top "slice_count" for 3D field, then majority thresholding
cfg.lamellarAngleDeg = 45;       %lamellar angle to x-axis
cfg.resultsRoot  = [];

fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn))
        cfg.(fn) = params.(fn);
    end
end

% Directories and folders

scriptDir   = fileparts(mfilename('fullpath'));
helpersDir  = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
customResultsRoot = ~isempty(cfg.resultsRoot);
if customResultsRoot
    resultsRoot = cfg.resultsRoot;
else
    resultsRoot = fullfile(scriptDir, 'results', 'sheets');
end
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

%-------------------------------

rng(cfg.rngSeed);
theta = deg2rad(cfg.lamellarAngleDeg);
coneBasis = [-cos(theta) sin(theta) 0; ...
             sin(theta)  cos(theta) 0; ...
             0           0          1];
[phi, meta] = spinodal_periodic_field(cfg.N, cfg.L, cfg.lambda_vox, cfg.bandwidth, cfg.nModes, cfg.coneDeg, coneBasis);

Svox = cfg.L / cfg.N;
tsV  = max(1, round(cfg.t_spin / Svox));
tbV  = max(3, round(cfg.t_base / Svox));
Nz   = tbV + tsV;

slab = min(max(1, round(cfg.slice_count)), size(phi,3));
tVol    = prctile(phi(:), 100*cfg.solid_frac);
fullMask = phi > tVol;
sliceMask = fullMask(:,:,max(1, end-slab+1):end);
if isempty(sliceMask)
    sliceMask = fullMask(:,:,end);
end
mask2 = mean(sliceMask, 3) >= 0.5;
% Force tileable boundaries in the base 2D mask so x/y edges coincide
mask2(end,:) = mask2(1,:);
mask2(:,end) = mask2(:,1);
solidFractionPattern = mean(mask2(:));

% Spin layer (tsV thick) atop the base; optionally etch away a fraction with
% a secondary GRF restricted to the spin region only.
spinLayer = repmat(mask2, 1, 1, tsV);
spinVoxelsInit = nnz(spinLayer);
removedFracActual = 0;
if spinVoxelsInit > 0 && cfg.remove_top_spin_frac > 0
    % Independent periodic GRF over the spin layer volume; isotropic cones.
    Nvec_rem = [cfg.N cfg.N tsV];
    Lvec_rem = [cfg.L cfg.L cfg.t_spin];
    rng(cfg.rngSeed + 17); % deterministic but separate stream
    phi_rem = spinodal_periodic_field_rect(Nvec_rem, Lvec_rem, ...
        cfg.lambda_vox, cfg.bandwidth, cfg.nModes, [90 90 90]);
    vals = phi_rem(spinLayer);
    cutoff = prctile(vals, 100*cfg.remove_top_spin_frac);
    removalMask = false(size(spinLayer));
    removalMask(spinLayer) = phi_rem(spinLayer) <= cutoff;
    % Enforce tileable boundaries on the removal mask to keep periodicity.
    removalMask(end,:,:) = removalMask(1,:,:);
    removalMask(:,end,:) = removalMask(:,1,:);
    spinLayer(removalMask) = false;
    % After removal, force spinLayer boundaries to match.
    spinLayer(end,:,:) = spinLayer(1,:,:);
    spinLayer(:,end,:) = spinLayer(:,1,:);
    removedFracActual = (spinVoxelsInit - nnz(spinLayer)) / spinVoxelsInit;
end

baseLayer = true(cfg.N, cfg.N, tbV);

sheetMask = cat(3, baseLayer, spinLayer);

tx = cfg.tilesXY(1); ty = cfg.tilesXY(2);
sheetMask = repmat(sheetMask, tx, ty, 1);

% Keep one fixed axis convention for all exports/pipeline steps.
sheetMask = permute(sheetMask, [2 1 3]);
Lx = cfg.L * ty;
Ly = cfg.L * tx;

solidFractionBase      = mean(baseLayer(:));
solidFractionSpinLayer = mean(spinLayer(:));
solidFractionSheet     = mean(sheetMask(:));
Lz = cfg.t_base + cfg.t_spin;

% Sanity check: edges must match for seamless tiling
edgeMismatchX = nnz(sheetMask(1,:,:) ~= sheetMask(end,:,:));
edgeMismatchY = nnz(sheetMask(:,1,:) ~= sheetMask(:,end,:));
if edgeMismatchX || edgeMismatchY
    error('PSSCone:periodicity', ...
        'Periodic tiling check failed (x mismatches: %d, y mismatches: %d).', ...
        edgeMismatchX, edgeMismatchY);
else
    fprintf('Periodicity check passed (x/y edges match).\n');
end

spinodalType = 'anisotropic';
if all(cfg.coneDeg >= 89)
    spinodalType = 'isotropic';
elseif cfg.coneDeg(1) == 30 && all(cfg.coneDeg([2 3]) == 0)
    spinodalType = 'lamellar';
elseif cfg.coneDeg(2) == 30 && all(cfg.coneDeg([1 3]) == 0)
    spinodalType = 'columnar';
elseif all(cfg.coneDeg == 15)
    spinodalType = 'Cubic';
end

runTimestamp = datetime('now');
tRatio = cfg.t_spin / max(cfg.t_base, eps);
ratioLabel = sprintf('tr%02d', round(100 * tRatio));
angleLabel = sprintf('ang%03d', round(cfg.lamellarAngleDeg));
runLabelBase = sprintf('sheetCone_%s_%s_%s_N%d_%dx%d', ...
    ratioLabel, angleLabel, lower(spinodalType), cfg.N, tx, ty);
if customResultsRoot
    typeRoot = resultsRoot;
else
    typeRoot = fullfile(resultsRoot, lower(spinodalType));
end
if ~exist(typeRoot, 'dir'), mkdir(typeRoot); end
runDir = unique_run_dir(typeRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

dx = cfg.L / cfg.N;
dz_base = cfg.t_base / tbV;
dz_spin = cfg.t_spin / tsV;

matPath = fullfile(runDir, 'sheet.mat');
[~, matBase] = fileparts(matPath);
voxelSizeXY = dx;
zVoxelThickness = [repmat(dz_base, 1, tbV), repmat(dz_spin, 1, tsV)];
save(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness', '-v7.3');
fprintf('Saved voxel mask to: %s\n', matPath);

manifest = struct();
manifest.mask = [matBase '.mat'];
manifest.var = 'sheetMask';
manifest.spacing = voxelSizeXY;
manifest.origin = [0 0 0];
manifest.material = 'SPINODAL';
manifest.notes = sprintf('Lamellar angle %.1f deg', cfg.lamellarAngleDeg);
manifestPath = fullfile(runDir, 'mesh_manifest.json');
try
    fidMan = fopen(manifestPath, 'w');
    if fidMan < 0
        warning('Unable to write manifest to %s', manifestPath);
    else
        fprintf(fidMan, '%s', jsonencode(manifest));
        fclose(fidMan);
        fprintf('Wrote manifest to: %s\n', manifestPath);
    end
catch ME
    warning('Failed to write manifest: %s', ME.message);
end

logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');
elapsedSeconds = toc(tStart);

paramLines = {
    sprintf('  grid: %d', cfg.N)
    sprintf('  cell_length_m: %.6f', cfg.L)
    sprintf('  lambda_vox: %.3f', cfg.lambda_vox)
    sprintf('  bandwidth: %.3f', cfg.bandwidth)
    sprintf('  nModes: %d', cfg.nModes)
    sprintf('  solid_frac_target: %.3f', cfg.solid_frac)
    sprintf('  coneDeg: [%s]', num2str(cfg.coneDeg))
    sprintf('  slice_count: %d', cfg.slice_count)
    sprintf('  lamellar_angle_deg: %.1f', cfg.lamellarAngleDeg)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    };

sheetLines = {
    sprintf('  tiles_xy: %dx%d', tx, ty)
    sprintf('  t_spin_mm: %.4f', cfg.t_spin*1e3)
    sprintf('  t_base_mm: %.4f', cfg.t_base*1e3)
    sprintf('  tsV/tbV/Nz (vox): %d / %d / %d', tsV, tbV, Nz)
    sprintf('  dims_mm: [%.3f %.3f %.3f]', Lx*1e3, Ly*1e3, Lz*1e3)
    sprintf('  pattern_frac_2d: %.3f', solidFractionPattern)
    sprintf('  solid_frac spin/base/sheet: %.3f / %.3f / %.3f', ...
            solidFractionSpinLayer, solidFractionBase, solidFractionSheet)
    sprintf('  spin_top_removed_target: %.3f', cfg.remove_top_spin_frac)
    sprintf('  spin_top_removed_actual: %.3f', removedFracActual)
    };

logLines = [
    {
    'Spinodoid sheet (rotated lamella) run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Folder: %s', runFolderName)
    sprintf('Type: %s (periodic XY)', spinodalType)
    sprintf('Mask MAT: %s.mat', matBase)
    ''
    'Base cell parameters:'
    } ;
    paramLines ;
    {
    ''
    'Sheet metadata:'
    } ;
    sheetLines
    {
    ''
    'Runtime:'
    sprintf('  elapsed_seconds: %.2f', elapsedSeconds)
    sprintf('  elapsed_minutes: %.2f', elapsedSeconds/60)
    }
    ];

fid = fopen(logPath, 'w');
if fid < 0
    warning('Unable to write log to %s', logPath);
else
    for i = 1:numel(logLines)
        fprintf(fid, '%s\n', logLines{i});
    end
    fclose(fid);
    fprintf('Logged run parameters to: %s\n', logPath);
end

fprintf('\n---- ROTATED SHEET SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Type: %s (sheet)\n', spinodalType);
fprintf('Grid: %d^3, lambda_vox: %.3f, bandwidth: %.3f\n', cfg.N, cfg.lambda_vox, cfg.bandwidth);
fprintf('Lamellar angle: %.1f deg\n', cfg.lamellarAngleDeg);
fprintf('Solid fraction (target / pattern / sheet): %.3f / %.3f / %.3f\n', ...
        cfg.solid_frac, solidFractionPattern, solidFractionSheet);
fprintf('Tiles: %dx%d, Mask MAT: %s\n', tx, ty, matPath);
fprintf('----------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
    fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end
