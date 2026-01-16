function PSSCone(params)
tStart = tic;
% PSSCone - Sheet generator with controllable in-plane lamellar angle.

if nargin < 1
    params = struct();
end

% ------------------------- Design knobs ---------------------------------------

cfg = struct();
cfg.N            = 64;          % grid size (NxNxN). Use powers of two for speed
cfg.L            = 40e-3;       % physical box length (mm)
cfg.lambda_vox   = 25;          % target feature wavelength in voxels (~rib/ligament spacing)
cfg.bandwidth    = 2;           % relative shell thickness around target |k| (0.1–0.3)
cfg.nModes       = 4000;        % number of Fourier modes to sample (1k–10k typical)
cfg.solid_frac   = 0.50;        % volume fraction of SOLID after threshold (0..1)
cfg.coneDeg      = [30 0 0];    % cone half-angles about x,y,z (90= isotropic). e.g. [90 90 90]
cfg.rngSeed      = 1;           % reproducible
cfg.sigma_vox    = 0.0;

cfg.t_spin       = 1e-3;        %spinodal thickness
cfg.t_base       = 2e-3;        %base thickness
cfg.tilesXY      = [1 1];       %tiling for periodicity
cfg.add_outer_skin_vox = 0;
cfg.slice_count  = 8;
cfg.align_with_cube = true;
cfg.lamellarAngleDeg = 0;       %lamellar angle to x-axis
cfg.resultsRoot  = [];

fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn))
        cfg.(fn) = params.(fn);
    end
end

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

rng(cfg.rngSeed);
[phi, meta] = spinodal_periodic_field(cfg.N, cfg.L, cfg.lambda_vox, cfg.bandwidth, cfg.nModes, cfg.coneDeg);

if cfg.sigma_vox > 0
    phi = periodic_gaussian_blur(phi, cfg.sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

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
mask2 = rotate_periodic_mask(mask2, cfg.lamellarAngleDeg);
solidFractionPattern = mean(mask2(:));

spinLayer = repmat(mask2, 1, 1, tsV);
baseLayer = true(cfg.N, cfg.N, tbV);

sheetMask = cat(3, baseLayer, spinLayer);

tx = cfg.tilesXY(1); ty = cfg.tilesXY(2);
sheetMask = repmat(sheetMask, tx, ty, 1);

if cfg.add_outer_skin_vox > 0
    s = min(cfg.add_outer_skin_vox, floor(size(sheetMask,1)/8));
    sheetMask(1:s,:,:)                 = true;
    sheetMask(end-s+1:end,:,:)         = true;
    sheetMask(:,1:s,:)                 = true;
    sheetMask(:,end-s+1:end,:)         = true;
end

if cfg.align_with_cube
    sheetMask = permute(sheetMask, [2 1 3]);
    Lx = cfg.L * ty;
    Ly = cfg.L * tx;
else
    Lx = cfg.L * tx;
    Ly = cfg.L * ty;
end

solidFractionBase      = mean(baseLayer(:));
solidFractionSpinLayer = mean(spinLayer(:));
solidFractionSheet     = mean(sheetMask(:));
Lz = cfg.t_base + cfg.t_spin;

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

stlPath = unique_path(fullfile(runDir, 'sheet.stl'));
[~, stlFileName, ext] = fileparts(stlPath);
stlFileName = [stlFileName ext];

sheetMaskPad = cat(1, false(1,size(sheetMask,2),size(sheetMask,3)), sheetMask, false(1,size(sheetMask,2),size(sheetMask,3)));
sheetMaskPad = cat(2, false(size(sheetMaskPad,1),1,size(sheetMaskPad,3)), sheetMaskPad, false(size(sheetMaskPad,1),1,size(sheetMaskPad,3)));
sheetMaskPad = cat(3, false(size(sheetMaskPad,1),size(sheetMaskPad,2),1), sheetMaskPad, false(size(sheetMaskPad,1),size(sheetMaskPad,2),1));

NxTot = size(sheetMask,1);
NyTot = size(sheetMask,2);
dx = cfg.L / cfg.N;
dy = dx;

dz_base = cfg.t_base / tbV;
dz_spin = cfg.t_spin / tsV;
zCore = [ (0:tbV-1)*dz_base + dz_base/2 , ...
          tbV*dz_base + (0:tsV-1)*dz_spin + dz_spin/2 ];

xCenters = ((0:NxTot-1) + 0.5) * dx;
yCenters = ((0:NyTot-1) + 0.5) * dy;

xCentersPad = [xCenters(1)-dx, xCenters, xCenters(end)+dx];
yCentersPad = [yCenters(1)-dy, yCenters, yCenters(end)+dy];
zCentersPad = [zCore(1)-dz_base, zCore, zCore(end)+dz_spin];

[X,Y,Z] = ndgrid(xCentersPad, yCentersPad, zCentersPad);
fv = isosurface(X, Y, Z, double(sheetMaskPad), 0.5);
if isempty(fv.vertices)
    error('Empty surface: adjust parameters to avoid trivial mask.');
end

meshStats.numFaces = size(fv.faces,1);
meshStats.numVertices = size(fv.vertices,1);

stlwrite_ascii(stlPath, fv.vertices, double(fv.faces));
fprintf('Wrote STL: %s\n', stlPath);

[matDir, matBase] = fileparts(stlPath);
matPath = fullfile(matDir, sprintf('%s.mat', matBase));
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
    sprintf('  sigma_vox: %.3f', cfg.sigma_vox)
    sprintf('  slice_count: %d', cfg.slice_count)
    sprintf('  lamellar_angle_deg: %.1f', cfg.lamellarAngleDeg)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    };

sheetLines = {
    sprintf('  tiles_xy: %dx%d', tx, ty)
    sprintf('  t_spin_mm: %.4f', cfg.t_spin*1e3)
    sprintf('  t_base_mm: %.4f', cfg.t_base*1e3)
    sprintf('  tsV/tbV/Nz (vox): %d / %d / %d', tsV, tbV, Nz)
    sprintf('  outer_skin_vox: %d', cfg.add_outer_skin_vox)
    sprintf('  dims_mm: [%.3f %.3f %.3f]', Lx*1e3, Ly*1e3, Lz*1e3)
    sprintf('  pattern_frac: %.3f', solidFractionPattern)
    sprintf('  solid_frac spin/base/sheet: %.3f / %.3f / %.3f', ...
            solidFractionSpinLayer, solidFractionBase, solidFractionSheet)
    };

logLines = [
    {
    'Spinodoid sheet (rotated lamella) run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Folder: %s', runFolderName)
    sprintf('Type: %s (periodic XY)', spinodalType)
    sprintf('STL: %s', stlFileName)
    sprintf('Faces: %d', meshStats.numFaces)
    sprintf('Vertices: %d', meshStats.numVertices)
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
fprintf('Tiles: %dx%d, STL: %s\n', tx, ty, stlPath);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('----------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
    fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end

function rotMask = rotate_periodic_mask(mask, angleDeg)
if abs(angleDeg) < 1e-6
    rotMask = mask;
    return;
end
N = size(mask,1);
[X,Y] = ndgrid(0:N-1, 0:N-1);
center = (N-1)/2;
Xc = X - center;
Yc = Y - center;
theta = deg2rad(angleDeg);
Xr = Xc*cos(theta) - Yc*sin(theta);
Yr = Xc*sin(theta) + Yc*cos(theta);
Xr = mod(Xr + center, N);
Yr = mod(Yr + center, N);
F = griddedInterpolant({0:N-1,0:N-1}, double(mask), 'nearest', 'nearest');
rotMask = F(Xr, Yr) >= 0.5;
end
