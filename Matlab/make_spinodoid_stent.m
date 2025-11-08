function make_spinodoid_stent(params)
%MAKE_SPINODOID_STENT Generate a periodic spinodoid stent wrapped on a cylinder.
%   make_spinodoid_stent() uses defaults; pass a struct to override fields.
%
%   Key parameters (defaults in parentheses):
%     Ri, Ro        : inner/outer radius [m] 
%     H             : axial span [m] 
%     Nt, Nz, Nr    : samples along theta, z, radial axes (256, 128, 48)
%     lambda_vox    : target wavelength in voxels (32)
%     bandwidth     : spectral shell width (0.25)
%     nModes        : number of Fourier modes (4000)
%     solid_frac    : target solid fraction (0.5)
%     coneDeg       : cone half-angles ( [90 90 90] )
%     rngSeed       : random seed (1)
%     sigma_vox     : Gaussian blur strength before threshold (0)
%     keepLargest   : keep largest component (true)
%     outer_skin_vox: optional axial skin after wrapping (0)

if nargin < 1
    params = struct();
end

cfg = struct( ...
    'Ri', 5e-3, ...
    'Ro', 6e-3, ...
    'H', 20e-3, ...
    'Nt', 256, ...
    'Nz', 128, ...
    'Nr', 48, ...
    'lambda_vox', 32, ...
    'bandwidth', 0.25, ...
    'nModes', 4000, ...
    'solid_frac', 0.50, ...
    'coneDeg', [90 90 90], ...
    'rngSeed', 1, ...
    'sigma_vox', 0.0, ...
    'keepLargest', true, ...
    'outer_skin_vox', 0);

fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn))
        cfg.(fn) = params.(fn);
    end
end

cfg.Nt = max(8, round(cfg.Nt));
cfg.Nz = max(4, round(cfg.Nz));
cfg.Nr = max(3, round(cfg.Nr));

if cfg.Ro <= cfg.Ri
    error('make_spinodoid_stent:invalidRadius', 'Outer radius must exceed inner radius.');
end

scriptDir  = fileparts(mfilename('fullpath'));
helpersDir = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir, 'results', 'stents');
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end
addpath(fullfile(scriptDir,'helpers'));

% ---- Generate periodic field in (theta,z,r) space -----------------------
rng(cfg.rngSeed);
Ltheta = 2*pi * ((cfg.Ri + cfg.Ro)/2);
Lz = cfg.H;
Lr = cfg.Ro - cfg.Ri;

[phi, meta] = spinodal_periodic_field_rect([cfg.Nt cfg.Nz cfg.Nr], ...
                                           [Ltheta Lz Lr], ...
                                           cfg.lambda_vox, cfg.bandwidth, ...
                                           cfg.nModes, cfg.coneDeg);

if cfg.sigma_vox > 0
    phi = periodic_gaussian_blur(phi, cfg.sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

t = prctile(phi(:), 100*cfg.solid_frac);
solidMask = phi > t;

if cfg.outer_skin_vox > 0
    s = min(cfg.outer_skin_vox, floor(cfg.Nz/8));
    solidMask(:,1:s,:)                 = true;
    solidMask(:,end-s+1:end,:)         = true;
end

axial_skin_vox = max(1, round(0.02*cfg.Nz));
axial_skin_vox = min(axial_skin_vox, max(1, floor(cfg.Nz/2)));
Nt = size(solidMask,1);
Nr = size(solidMask,3);
for i = 1:Nt
    for r = 1:Nr
        column = squeeze(solidMask(i,:,r));
        firstSolid = find(column, 1, 'first');
        if ~isempty(firstSolid) && firstSolid <= axial_skin_vox + 1
            column(1:firstSolid) = true;
        end
        lastSolid = find(column, 1, 'last');
        if ~isempty(lastSolid) && lastSolid >= (cfg.Nz - axial_skin_vox)
            column(lastSolid:end) = true;
        end
        solidMask(i,:,r) = column;
    end
end

if cfg.keepLargest && exist('bwconncomp','file') == 2
    try
        maskWrap = cat(1, solidMask, solidMask);
        CC = bwconncomp(maskWrap, 26);
        [~,iMax] = max(cellfun(@numel, CC.PixelIdxList));
        tmp = false(size(maskWrap));
        tmp(CC.PixelIdxList{iMax}) = true;
        solidMask = tmp(1:size(solidMask,1),:,:);
    catch
        % if toolbox missing, keep original mask
    end
end

solidMask(:,end,:) = solidMask(:,1,:); % ensure periodicity in theta

solidFractionActual = mean(solidMask(:));

% ---- Build run folder ---------------------------------------------------
spinodalType = 'anisotropic';
if all(cfg.coneDeg >= 89)
    spinodalType = 'isotropic';
end

runTimestamp = datetime('now');
runLabelBase = sprintf('stent_%s_Ri%03dum_Ro%03dum_H%03dum', ...
    lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), round(1e6*cfg.H));
runDir = unique_run_dir(resultsRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('spinodoid_STENT_Ri%03dum_Ro%03dum_H%03dum.stl', ...
    round(1e6*cfg.Ri), round(1e6*cfg.Ro), round(1e6*cfg.H));
stlPath = unique_path(fullfile(runDir, stlFileName));

% ---- Extract iso-surface in physical space (stent-copy style) -----------
thetaVec = linspace(0, 2*pi, cfg.Nt+1); thetaVec(end) = [];
if cfg.Nz > 1
    zVec = linspace(0, cfg.H, cfg.Nz);
else
    zVec = 0;
end
if cfg.Nr > 1
    rVec = linspace(cfg.Ri, cfg.Ro, cfg.Nr);
else
    rVec = cfg.Ri;
end

maskRTZ = permute(solidMask, [3 1 2]); % R x theta x z
[maskRadial, rPad] = padRadial(maskRTZ, rVec);
[maskAxial, zPad] = padAxial(maskRadial, zVec);
[maskWrapped, thetaWrap] = wrapTheta(maskAxial, thetaVec);
[XGrid, YGrid, ZGrid] = cylindricalGrid(rPad, thetaWrap, zPad);

fv = isosurface(XGrid, YGrid, ZGrid, double(maskWrapped), 0.5);
if isempty(fv.vertices)
    error('make_spinodoid_stent:emptySurface', 'Stent surface is empty. Check parameters.');
end

fv.faces = double(fv.faces);

% Ensure face orientation matches outward normal (right-hand rule)
N = faceNormals(fv.vertices, fv.faces);
centroids = (fv.vertices(fv.faces(:,1),:) + fv.vertices(fv.faces(:,2),:) + fv.vertices(fv.faces(:,3),:))/3;
radial = centroids(:,1).^2 + centroids(:,2).^2;
flipIdx = radial < ((cfg.Ri + cfg.Ro)/2)^2 & dot(N, centroids, 2) > 0;
fv.faces(flipIdx,:) = fv.faces(flipIdx,[1 3 2]);

meshStats.numFaces = size(fv.faces,1);
meshStats.numVertices = size(fv.vertices,1);

stlwrite_ascii(stlPath, fv.vertices, fv.faces);
fprintf('Wrote STL: %s\n', stlPath);

% ---- Logging ------------------------------------------------------------
logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');

paramLines = {
    sprintf('  Ri_mm: %.3f', cfg.Ri*1e3)
    sprintf('  Ro_mm: %.3f', cfg.Ro*1e3)
    sprintf('  H_mm: %.3f', cfg.H*1e3)
    sprintf('  Nt/Nz/Nr: %d / %d / %d', cfg.Nt, cfg.Nz, cfg.Nr)
    sprintf('  lambda_vox: %.3f', cfg.lambda_vox)
    sprintf('  bandwidth: %.3f', cfg.bandwidth)
    sprintf('  nModes: %d', cfg.nModes)
    sprintf('  solid_frac_target: %.3f', cfg.solid_frac)
    sprintf('  coneDeg: [%s]', num2str(cfg.coneDeg))
    sprintf('  sigma_vox: %.3f', cfg.sigma_vox)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    };

geomLines = {
    sprintf('  shell_thickness_mm: %.3f', (cfg.Ro-cfg.Ri)*1e3)
    sprintf('  outer_skin_vox: %d', cfg.outer_skin_vox)
    sprintf('  solid_frac_actual: %.3f', solidFractionActual)
    sprintf('  modes_used: %d', meta.nModes)
    };

logLines = [
    {
    'Spinodoid stent run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Folder: %s', runFolderName)
    sprintf('Type: %s', spinodalType)
    sprintf('STL: %s', stlFileName)
    sprintf('Faces: %d', meshStats.numFaces)
    sprintf('Vertices: %d', meshStats.numVertices)
    ''
    'Parameters:'
    } ;
    paramLines ;
    {
    ''
    'Geometry / stats:'
    } ;
    geomLines
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

fprintf('\n---- STENT RUN SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Ri/Ro/H (mm): %.3f / %.3f / %.3f\n', cfg.Ri*1e3, cfg.Ro*1e3, cfg.H*1e3);
fprintf('Grid Nt/Nz/Nr: %d / %d / %d\n', cfg.Nt, cfg.Nz, cfg.Nr);
fprintf('Solid fraction (target/actual): %.3f / %.3f\n', cfg.solid_frac, solidFractionActual);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('STL path: %s\n', stlPath);
fprintf('---------------------------\n\n');
end

function [maskPad, rPad] = padRadial(mask, rVec)
dr = spacing_stent(rVec);
rPad = [max(rVec(1) - dr, 0), rVec(:).', rVec(end) + dr];
maskPad = false(size(mask) + [2 0 0]);
maskPad(2:end-1,:,:) = mask;
end

function [maskPad, zPad] = padAxial(mask, zVec)
dz = spacing_stent(zVec);
zPad = [zVec(1) - dz, zVec(:).', zVec(end) + dz];
maskPad = false(size(mask) + [0 0 2]);
maskPad(:,:,2:end-1) = mask;
end

function [maskWrap, thetaWrap] = wrapTheta(mask, thetaVec)
dtheta = spacing_stent(thetaVec);
thetaWrap = [thetaVec(:).', thetaVec(end) + dtheta];
maskWrap = cat(2, mask, mask(:,1,:));
end

function [X,Y,Z] = cylindricalGrid(rVec, thetaVec, zVec)
[R,Theta,Z] = ndgrid(rVec, thetaVec, zVec);
X = R .* cos(Theta);
Y = R .* sin(Theta);
end

function d = spacing_stent(vec)
vec = vec(:).';
if numel(vec) > 1
    d = mean(diff(vec));
else
    d = max(vec,1);
end
end
