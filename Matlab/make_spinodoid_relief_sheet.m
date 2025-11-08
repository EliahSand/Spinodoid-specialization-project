function make_spinodoid_relief_sheet(params)
%MAKE_SPINODOID_RELIEF_SHEET Solid slab + extruded spinodoid relief.
%   Bottom: dense homogeneous slab (thickness t_base)
%   Top   : constant-thickness spinodoid "carpet" extruded from a 2D slice
%           of the same periodic field used in main.m (no taper/slope).
%   The two layers are fused into a single watertight solid block.

if nargin < 1
    params = struct();
end

%% -------------------- Design knobs (with defaults) ---------------------
cfg.N            = 128;        % grid size per base cell
cfg.L            = 10e-3;      % physical cell size [m]
cfg.lambda_vox   = 50;
cfg.bandwidth    = 0.22;
cfg.nModes       = 4000;
cfg.solid_frac   = 0.50;
cfg.coneDeg      = [90 90 90];
cfg.rngSeed      = 1;
cfg.sigma_vox    = 0.0;        % optional Gaussian blur strength

cfg.t_spin       = 2.0e-3;     % spinodoid relief thickness [m]
cfg.t_base       = 1.5e-3;     % dense base thickness [m]
cfg.tilesXY      = [1 1];      % tiling in X/Y
cfg.add_outer_skin_vox = 0;    % perimeter sealing on the whole sheet
cfg.slice_count  = 4;          % number of top slices to average for 2D driver

fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn)
        cfg.(fn) = params.(fn);
    end
end

%% -------------------- Paths --------------------------------------------
scriptDir   = fileparts(mfilename('fullpath'));
helpersDir  = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir, 'results', 'relief_sheets');
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

%% -------------------- Generate periodic spinodoid cell -----------------
rng(cfg.rngSeed);
[phi, meta] = spinodal_periodic_field(cfg.N, cfg.L, cfg.lambda_vox, ...
                                      cfg.bandwidth, cfg.nModes, cfg.coneDeg);

if cfg.sigma_vox > 0
    phi = periodic_gaussian_blur(phi, cfg.sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

%% -------------------- 2D driver & extrusion ----------------------------
slab = min(max(1, round(cfg.slice_count)), size(phi,3));
phi2 = mean(phi(:,:,1:slab), 3);           % average a thin top slab for stability
phi2 = phi2 - mean(phi2(:));
phi2 = phi2 ./ (std(phi2(:)) + eps);

t2      = prctile(phi2(:), 100*cfg.solid_frac);
mask2   = phi2 > t2;                        % 2D spinodal islands (periodic in XY)
solidFractionCell = mean(mask2(:));         % report 2D coverage

%% -------------------- Build layered sheet ------------------------------
Svox = cfg.L / cfg.N;
tsV  = max(1, round(cfg.t_spin / Svox));
tbV  = max(3, round(cfg.t_base / Svox));   % ensure plate survives marching cubes
Nz   = tbV + tsV;

spinLayer = repmat(mask2, 1, 1, tsV);       % pure extrusion: vertical walls, closed top
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

solidFractionSpinLayer = mean(spinLayer(:));
solidFractionBase      = mean(baseLayer(:));
solidFractionSheet     = mean(sheetMask(:));

Lx = cfg.L * tx;
Ly = cfg.L * ty;
Lz = cfg.t_base + cfg.t_spin;

%% -------------------- Export & log ------------------------------------
spinodalType = 'anisotropic';
if all(cfg.coneDeg >= 89)
    spinodalType = 'isotropic';
end

runTimestamp = datetime('now');
runLabelBase = sprintf('sheetRelief_%s_N%d_%dx%d_t%dum_%dum', ...
    lower(spinodalType), cfg.N, tx, ty, round(1e6*cfg.t_spin), round(1e6*cfg.t_base));
runDir = unique_run_dir(resultsRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('sheet_relief_spinodoid_N%d_%dx%d_L%.1fmm_t%.2f_%.2fmm.stl', ...
    cfg.N, tx, ty, cfg.L*1e3, cfg.t_spin*1e3, cfg.t_base*1e3);
stlPath = unique_path(fullfile(runDir, stlFileName));

% Pad with one-voxel void on all sides to ensure watertight isosurface
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

logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');

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
    sprintf('  rng_seed: %d', cfg.rngSeed)
    };

sheetLines = {
    sprintf('  tiles_xy: %dx%d', tx, ty)
    sprintf('  t_spin_mm: %.4f', cfg.t_spin*1e3)
    sprintf('  t_base_mm: %.4f', cfg.t_base*1e3)
    sprintf('  tsV/tbV/Nz (vox): %d / %d / %d', tsV, tbV, Nz)
    sprintf('  outer_skin_vox: %d', cfg.add_outer_skin_vox)
    sprintf('  dims_mm: [%.3f %.3f %.3f]', Lx*1e3, Ly*1e3, Lz*1e3)
    sprintf('  solid_frac spin/base/sheet: %.3f / %.3f / %.3f', ...
            solidFractionSpinLayer, solidFractionBase, solidFractionSheet)
    };

logLines = [
    {
    'Spinodoid relief sheet run log'
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
    ];

fid = fopen(logPath, 'w');
if fid < 0
    warning('Unable to write run log to %s', logPath);
else
    for i = 1:numel(logLines)
        fprintf(fid, '%s\n', logLines{i});
    end
    fclose(fid);
    fprintf('Logged run parameters to: %s\n', logPath);
end

fprintf('\n---- RUN SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Generated: %s\n', timestampStr);
fprintf('Type: %s (extruded relief sheet)\n', spinodalType);
fprintf('Grid: %d^3, lambda_vox: %.3f, bandwidth: %.3f\n', cfg.N, cfg.lambda_vox, cfg.bandwidth);
fprintf('Solid fraction (target / cell / sheet): %.3f / %.3f / %.3f\n', ...
        cfg.solid_frac, solidFractionCell, solidFractionSheet);
fprintf('Tiles: %dx%d, STL: %s\n', tx, ty, stlPath);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('----------------------\n\n');
end
