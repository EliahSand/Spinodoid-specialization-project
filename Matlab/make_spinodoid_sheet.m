function make_spinodoid_sheet()
% Two-layer sheet that shares the same field-generation pipeline as main.m:
% - Generate a periodic spinodoid cell via cosine-mode synthesis
% - Use the top tsV slices as the porous layer
% - Stack on top of a dense base slab, tile in X/Y, export STL + log

%% -------------------- Design knobs (main-style core) -------------------
N            = 128;        % grid size for the periodic cell (NxNxN)
L            = 40e-3;      % physical side length of one cell [m]
lambda_vox   = 32;         % target wavelength [vox]
bandwidth    = 0.22;       % spectral shell thickness
nModes       = 4000;       % number of Fourier modes
solid_frac   = 0.50;       % SOLID fraction (controls island coverage)
coneDeg      = [90 90 90]; % isotropic relief
rngSeed      = 1;          % reproducible field

% Optional smoothing prior to threshold
sigma_vox    = 0.8;        % Gaussian blur strength (0 disables)

% Sheet-specific geometry
tilesXY      = [2 2];      % replicate cell in X,Y
t_spin       = 2.0e-3;     % spinodoid layer thickness [m]
t_base       = 1.5e-3;     % dense base thickness [m]
add_outer_skin_vox = 0;    % perimeter sealing on full sheet (0 disables)

%% -------------------- Paths --------------------------------------------
scriptDir   = fileparts(mfilename('fullpath'));
helpersDir  = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir, 'results', 'sheets');
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

%% -------------------- Generate periodic field --------------------------
rng(rngSeed);
[phi, meta] = spinodal_periodic_field(N, L, lambda_vox, bandwidth, nModes, coneDeg);

if sigma_vox > 0
    phi = periodic_gaussian_blur(phi, sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

t         = prctile(phi(:), 100*solid_frac);
solidMask = phi > t;

solidFractionCell = mean(solidMask(:));

%% -------------------- Build layered sheet ------------------------------
Svox = L / N;                                % voxel size [m]
tsV  = max(1, round(t_spin / Svox));
tbV  = max(1, round(t_base / Svox));
Nz   = tbV + tsV;

if tsV > size(solidMask,3)
    error('Spin layer thickness exceeds available cell slices. Increase N or reduce t_spin.');
end

spinLayer = solidMask(:,:,1:tsV);
baseLayer = true(N, N, tbV);

sheetMask = cat(3, baseLayer, spinLayer);

tx = tilesXY(1); ty = tilesXY(2);
sheetMask = repmat(sheetMask, tx, ty, 1);

if add_outer_skin_vox > 0
    s = min(add_outer_skin_vox, floor(size(sheetMask,1)/8));
    sheetMask(1:s,:,:)                 = true;
    sheetMask(end-s+1:end,:,:)         = true;
    sheetMask(:,1:s,:)                 = true;
    sheetMask(:,end-s+1:end,:)         = true;
end

solidFractionSpinLayer = mean(spinLayer(:));
solidFractionBase      = mean(baseLayer(:));
solidFractionSheet     = mean(sheetMask(:));

Lx = L * tx;
Ly = L * ty;
Lz = Nz * Svox;

%% -------------------- Export & log ------------------------------------
spinodalType = 'anisotropic';
if all(coneDeg >= 89)
    spinodalType = 'isotropic';
end

runTimestamp = datetime('now');
runLabelBase = sprintf('sheet_%s_N%d_%dx%d_t%dum_%dum', ...
    lower(spinodalType), N, tx, ty, round(1e6 * t_spin), round(1e6 * t_base));
runDir = unique_run_dir(resultsRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('sheet_spinodoid_base_N%d_%dx%d_L%.1fmm_t%.2f_%.2fmm.stl', ...
    N, tx, ty, L*1e3, t_spin*1e3, t_base*1e3);
stlPath = unique_path(fullfile(runDir, stlFileName));

meshStats = binary_to_watertight_stl(sheetMask, Lx, stlPath);
fprintf('Wrote STL: %s\n', stlPath);

logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');

paramLines = {
    sprintf('  grid: %d', N)
    sprintf('  cell_length_m: %.6f', L)
    sprintf('  lambda_vox: %.3f', lambda_vox)
    sprintf('  bandwidth: %.3f', bandwidth)
    sprintf('  nModes: %d', nModes)
    sprintf('  solid_frac_target: %.3f', solid_frac)
    sprintf('  coneDeg: [%s]', num2str(coneDeg))
    sprintf('  sigma_vox: %.3f', sigma_vox)
    sprintf('  rng_seed: %d', rngSeed)
    };

sheetLines = {
    sprintf('  tiles_xy: %dx%d', tx, ty)
    sprintf('  t_spin_mm: %.4f', t_spin*1e3)
    sprintf('  t_base_mm: %.4f', t_base*1e3)
    sprintf('  tsV/tbV/Nz (vox): %d / %d / %d', tsV, tbV, Nz)
    sprintf('  outer_skin_vox: %d', add_outer_skin_vox)
    sprintf('  dims_mm: [%.3f %.3f %.3f]', Lx*1e3, Ly*1e3, Lz*1e3)
    sprintf('  solid_frac spin/base/sheet: %.3f / %.3f / %.3f', ...
            solidFractionSpinLayer, solidFractionBase, solidFractionSheet)
    };

logLines = [
    {
    'Spinodoid sheet run log'
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
fprintf('Type: %s (sheet)\n', spinodalType);
fprintf('Grid: %d^3, lambda_vox: %.3f, bandwidth: %.3f\n', N, lambda_vox, bandwidth);
fprintf('Solid fraction (target / cell / sheet): %.3f / %.3f / %.3f\n', ...
        solid_frac, solidFractionCell, solidFractionSheet);
fprintf('Tiles: %dx%d, STL: %s\n', tx, ty, stlPath);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('----------------------\n\n');
end
