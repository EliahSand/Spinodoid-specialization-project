%MAKE_LAMINATE_SHEET Generate a two-layer laminate mask and exports.
%   Bottom layer is fully solid.
%   Top layer is line-ribbed in X or Y.
%
%   Design rule enforced in the generated geometry:
%     bottom thickness = top thickness = line width
%
%   Inputs are interpreted in SI units (meters):
%     L [m], a [m]


function make_laminate_sheet(params)

tStart = tic;

if nargin < 1
    params = struct();
end

cfg = struct();
cfg.N = 128;                      % Voxels per side in one tile (XY)
cfg.L = 40e-3;                    % Physical cell length in X/Y for one tile (m)
cfg.a = 1e-3;                     % Target feature size (m): tb = ts = line width = a
cfg.lineDirection = 'y';          % 'x' (lines run along X) or 'y' (along Y)
cfg.lineGapFactor = 1.0;          % Gap width = lineGapFactor * line width
cfg.tilesXY = [10 1];              % Tile counts in X and Y
cfg.resultsRoot = [];             % Optional override

fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn))
        cfg.(fn) = params.(fn);
    end
end

cfg.lineDirection = validatestring(cfg.lineDirection, {'x', 'y'});
if cfg.N < 4
    error('N must be >= 4.');
end

scriptDir = fileparts(mfilename('fullpath'));
matlabRoot = fileparts(scriptDir);
helpersDir = fullfile(matlabRoot, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);

if isempty(cfg.resultsRoot)
    resultsRoot = fullfile(matlabRoot, 'results', 'sheets', 'laminate');
else
    resultsRoot = cfg.resultsRoot;
end
if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end

% Isotropic XY voxel size in meters.
Svox = cfg.L / cfg.N;

% Convert target size to an integer voxel count so equality is exact.
aV = max(1, round(cfg.a / Svox));
tbV = aV;
tsV = aV;
lineWv = aV;
gapV = max(1, round(cfg.lineGapFactor * lineWv));
pitchV = lineWv + gapV;

% Actual resolved geometry after voxelization.
aResolved = aV * Svox;
tBase = tbV * Svox;
tTop = tsV * Svox;

tx = cfg.tilesXY(1);
ty = cfg.tilesXY(2);
NxTot = cfg.N * tx;
NyTot = cfg.N * ty;

% Generate ribs over the full tiled extent to avoid seam thickening at tile joins.
if strcmpi(cfg.lineDirection, 'x')
    on = mod((0:NyTot-1), pitchV) < lineWv;
    lineMask2D = repmat(on, NxTot, 1);
else
    on = mod((0:NxTot-1), pitchV) < lineWv;
    lineMask2D = repmat(on(:), 1, NyTot);
end

baseLayer = true(NxTot, NyTot, tbV);
topLayer = repmat(lineMask2D, 1, 1, tsV);
sheetMask = cat(3, baseLayer, topLayer);

runTimestamp = datetime('now');
dirTag = lower(cfg.lineDirection);
aTagUm = round(1e6 * aResolved);
runLabelBase = sprintf('laminate_dir%s_a%04dum_N%d_%dx%d', dirTag, aTagUm, cfg.N, tx, ty);
runDir = unique_run_dir(resultsRoot, runLabelBase);

stlPath = unique_path(fullfile(runDir, 'sheet.stl'));
[~, stlBase] = fileparts(stlPath);
matPath = fullfile(runDir, [stlBase '.mat']);

% Build isosurface with anisotropic Z spacing support.
dzBase = tBase / tbV;
dzTop = tTop / tsV;
zVoxelThickness = [repmat(dzBase, 1, tbV), repmat(dzTop, 1, tsV)];

sheetMaskPad = cat(1, false(1, size(sheetMask,2), size(sheetMask,3)), sheetMask, false(1, size(sheetMask,2), size(sheetMask,3)));
sheetMaskPad = cat(2, false(size(sheetMaskPad,1), 1, size(sheetMaskPad,3)), sheetMaskPad, false(size(sheetMaskPad,1), 1, size(sheetMaskPad,3)));
sheetMaskPad = cat(3, false(size(sheetMaskPad,1), size(sheetMaskPad,2), 1), sheetMaskPad, false(size(sheetMaskPad,1), size(sheetMaskPad,2), 1));

dx = Svox;
dy = Svox;
zCore = [ (0:tbV-1)*dzBase + dzBase/2 , ...
          tbV*dzBase + (0:tsV-1)*dzTop + dzTop/2 ];

xCenters = ((0:NxTot-1) + 0.5) * dx;
yCenters = ((0:NyTot-1) + 0.5) * dy;
xCentersPad = [xCenters(1)-dx, xCenters, xCenters(end)+dx];
yCentersPad = [yCenters(1)-dy, yCenters, yCenters(end)+dy];
zCentersPad = [zCore(1)-dzBase, zCore, zCore(end)+dzTop];

[X, Y, Z] = ndgrid(xCentersPad, yCentersPad, zCentersPad);
fv = isosurface(X, Y, Z, double(sheetMaskPad), 0.5);
if isempty(fv.vertices)
    error('Empty laminate surface; adjust N or a.');
end
stlwrite_ascii(stlPath, fv.vertices, double(fv.faces));

voxelSizeXY = Svox;
save(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness', '-v7.3');

manifest = struct();
manifest.mask = [stlBase '.mat'];
manifest.var = 'sheetMask';
manifest.spacing = voxelSizeXY;
manifest.origin = [0 0 0];
manifest.material = 'SPINODAL';
manifest.notes = sprintf('Laminate lines along %s', upper(cfg.lineDirection));
manifestPath = fullfile(runDir, 'mesh_manifest.json');
fidMan = fopen(manifestPath, 'w');
if fidMan >= 0
    fprintf(fidMan, '%s', jsonencode(manifest));
    fclose(fidMan);
else
    warning('Could not write manifest: %s', manifestPath);
end

elapsedSeconds = toc(tStart);
logPath = fullfile(runDir, 'run_log.txt');
fid = fopen(logPath, 'w');
if fid >= 0
    fprintf(fid, '=== Laminate Sheet Generation Log ===\n');
    fprintf(fid, 'timestamp: %s\n', datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS'));
    fprintf(fid, 'run_dir: %s\n', runDir);
    fprintf(fid, 'N: %d\n', cfg.N);
    fprintf(fid, 'L_m: %.6e\n', cfg.L);
    fprintf(fid, 'voxel_size_xy_m: %.6e\n', Svox);
    fprintf(fid, 'line_direction: %s\n', cfg.lineDirection);
    fprintf(fid, 'line_gap_factor: %.3f\n', cfg.lineGapFactor);
    fprintf(fid, 'a_target_m: %.6e\n', cfg.a);
    fprintf(fid, 'a_resolved_m: %.6e\n', aResolved);
    fprintf(fid, 't_base_mm: %.6f\n', 1e3 * tBase);
    fprintf(fid, 't_spin_mm: %.6f\n', 1e3 * tTop);
    fprintf(fid, 't_base_m: %.6e\n', tBase);
    fprintf(fid, 't_spin_m: %.6e\n', tTop);
    fprintf(fid, 'line_width_m: %.6e\n', aResolved);
    fprintf(fid, 'line_width_vox: %d\n', lineWv);
    fprintf(fid, 'base_slices: %d\n', tbV);
    fprintf(fid, 'spin_slices: %d\n', tsV);
    fprintf(fid, 'tiles_xy: [%d %d]\n', tx, ty);
    fprintf(fid, 'elapsed_s: %.3f\n', elapsedSeconds);
    fclose(fid);
end

fprintf('Laminate model generated.\n');
fprintf('  Run folder: %s\n', runDir);
fprintf('  STL: %s\n', stlPath);
fprintf('  MAT: %s\n', matPath);
fprintf('  Geometry rule: t_base = t_top = line_width = %.6e m (%d vox)\n', aResolved, aV);

end
