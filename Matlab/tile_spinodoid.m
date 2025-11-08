function tile_spinodoid()
% Multi-cell spinodoid generator built on the same pipeline as main.m
% 1) Generate a single periodic field via cosine-mode synthesis
% 2) Threshold, optionally clean/skin the cell
% 3) Tile across X/Y/Z and export STL + run log

%% -------------------- Design knobs (main-style) ------------------------
N            = 128;        % grid size of the base cell (NxNxN)
L            = 0.50;       % physical box length of one cell (arbitrary units)
lambda_vox   = 32;         % target feature wavelength (in voxels)
bandwidth    = 0.30;       % spectral shell thickness (0.1–0.3)
nModes       = 4000;       % number of Fourier modes
solid_frac   = 0.50;       % SOLID volume fraction after threshold (0..1)
coneDeg      = [90 90 90]; % cone half-angles; [90 90 90] = isotropic
rngSeed      = 1;          % reproducible seed
skin_ratio   = 0.02;       % exterior skin thickness as fraction of N (per cell)
keepLargest  = true;       % keep only the largest connected solid component

% Optional smoothing of phi BEFORE threshold (still periodic)
sigma_vox    = 0.0;        % Gaussian blur (0 disables)

% Tiling of the periodic cell
tiles        = [3 3 1];    % [tx ty tz] or scalar shorthand
outer_skin_vox = 0;        % extra skin on the final tiled domain (0 disables)

%% -------------------- Derived paths ------------------------------------
scriptDir   = fileparts(mfilename('fullpath'));
helpersDir  = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir, 'results', 'tiled');
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

%% -------------------- Generate base periodic cell ----------------------
rng(rngSeed);
[phi, meta] = spinodal_periodic_field(N, L, lambda_vox, bandwidth, nModes, coneDeg);

if sigma_vox > 0
    phi = periodic_gaussian_blur(phi, sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

t         = prctile(phi(:), 100*solid_frac);
solidMask = phi > t;

skin_thickness_vox = max(0, round(skin_ratio * N));
if skin_thickness_vox > 0
    s = skin_thickness_vox;
    solidMask(1:s,:,:)                 = true;
    solidMask(end-s+1:end,:,:)         = true;
    solidMask(:,1:s,:)                 = true;
    solidMask(:,end-s+1:end,:)         = true;
    solidMask(:,:,1:s)                 = true;
    solidMask(:,:,end-s+1:end)         = true;
end

if keepLargest
    try
        if exist('bwconncomp','file') == 2
            CC = bwconncomp(solidMask, 26);
            [~,iMax] = max(cellfun(@numel, CC.PixelIdxList));
            tmp = false(size(solidMask));
            tmp(CC.PixelIdxList{iMax}) = true;
            solidMask = tmp;
        end
    catch
        % if toolbox missing, skip silently
    end
end

solidFractionCell = mean(solidMask(:));

%% -------------------- Tile the domain ----------------------------------
if isscalar(tiles), tiles = [tiles tiles tiles]; end
tx = tiles(1); ty = tiles(2); tz = tiles(3);

tiledMask = repmat(solidMask, tx, ty, tz);
Lbig = L * [tx ty tz];

if outer_skin_vox > 0
    s = min(outer_skin_vox, floor(size(tiledMask,1)/8));
    tiledMask(1:s,:,:)                 = true;
    tiledMask(end-s+1:end,:,:)         = true;
    tiledMask(:,1:s,:)                 = true;
    tiledMask(:,end-s+1:end,:)         = true;
    tiledMask(:,:,1:s)                 = true;
    tiledMask(:,:,end-s+1:end)         = true;
end

solidFractionTiled = mean(tiledMask(:));

%% -------------------- Export & log ------------------------------------
spinodalType = 'anisotropic';
if all(coneDeg >= 89)
    spinodalType = 'isotropic';
end

runTimestamp = datetime('now');
runLabelBase = sprintf('tiled_%s_N%d_sf%02d_%dx%dx%d', ...
    lower(spinodalType), N, round(100*solid_frac), tx, ty, tz);
runDir = unique_run_dir(resultsRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('spinodoid_TILED_N%d_sf%02d_%dx%dx%d.stl', ...
    N, round(100*solid_frac), tx, ty, tz);
stlPath = unique_path(fullfile(runDir, stlFileName));

meshStats = binary_to_watertight_stl(tiledMask, Lbig(1), stlPath);
fprintf('Wrote STL: %s\n', stlPath);

logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');

paramLines = {
    sprintf('  grid: %d', N)
    sprintf('  cell_length: %.5f', L)
    sprintf('  lambda_vox: %.3f', lambda_vox)
    sprintf('  bandwidth: %.3f', bandwidth)
    sprintf('  nModes: %d', nModes)
    sprintf('  solid_frac_target: %.3f', solid_frac)
    sprintf('  coneDeg: [%s]', num2str(coneDeg))
    sprintf('  skin_ratio: %.3f (vox=%d)', skin_ratio, skin_thickness_vox)
    sprintf('  keepLargest: %d', keepLargest)
    sprintf('  sigma_vox: %.3f', sigma_vox)
    sprintf('  rng_seed: %d', rngSeed)
    };

tileLines = {
    sprintf('  tiles: %dx%dx%d', tx, ty, tz)
    sprintf('  L_big: [%.5f %.5f %.5f]', Lbig)
    sprintf('  outer_skin_vox: %d', outer_skin_vox)
    sprintf('  solid_frac_cell/tiled: %.3f / %.3f', ...
            solidFractionCell, solidFractionTiled)
    };

logLines = [
    {
    'Spinodoid tiled run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Folder: %s', runFolderName)
    sprintf('Type: %s (periodic base cell)', spinodalType)
    sprintf('STL: %s', stlFileName)
    sprintf('Faces: %d', meshStats.numFaces)
    sprintf('Vertices: %d', meshStats.numVertices)
    ''
    'Base cell parameters:'
    } ;
    paramLines ;
    {
    ''
    'Tiling metadata:'
    } ;
    tileLines
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
fprintf('Type: %s\n', spinodalType);
fprintf('Grid: %d^3, lambda_vox: %.3f, bandwidth: %.3f\n', N, lambda_vox, bandwidth);
fprintf('Solid fraction (target / cell / tiled): %.3f / %.3f / %.3f\n', ...
        solid_frac, solidFractionCell, solidFractionTiled);
fprintf('Tiles: %dx%dx%d, STL: %s\n', tx, ty, tz, stlPath);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('----------------------\n\n');
end
