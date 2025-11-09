function tile_spinodoid_stent(params)
tStart = tic;
%TILE_SPINODOID_STENT Generate multiple spinodoid stents fused along the axis.

if nargin < 1, params = struct(); end

cfg = struct('Ri',5e-3,'Ro',6e-3,'H',20e-3,'Nt',256,'Nz',128,'Nr',48, ...
             'lambda_vox',32,'bandwidth',0.25,'nModes',4000,'solid_frac',0.5, ...
             'coneDeg',[90 90 90],'rngSeed',1,'sigma_vox',0.0,'keepLargest',true, ...
             'outer_skin_vox',0,'tilesAxial',3);
fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn)), cfg.(fn) = params.(fn); end
end

cfg.tilesAxial = max(1, round(cfg.tilesAxial));
cfg.Nt = max(8, round(cfg.Nt));
cfg.Nz = max(4, round(cfg.Nz));
cfg.Nr = max(3, round(cfg.Nr));
if cfg.Ro <= cfg.Ri, error('tile_spinodoid_stent:invalidRadius', 'Ro must exceed Ri'); end

scriptDir = fileparts(mfilename('fullpath'));
helpersDir = fullfile(scriptDir,'helpers');
if ~isfolder(helpersDir), error('Helpers directory missing: %s', helpersDir); end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir,'results','stents_tiled');
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

[solidMaskBase, meta] = generate_stent_mask(cfg);
solidMask = repmat(solidMaskBase, [1 cfg.tilesAxial 1]);
NzTotal = size(solidMask,2);
HTotal = cfg.H * cfg.tilesAxial;

% Optional outer (axial) skin at both ends of the tiled stack
if cfg.outer_skin_vox > 0
    s = min(cfg.outer_skin_vox, floor(NzTotal/8));
    solidMask(:,1:s,:) = true;
    solidMask(:,end-s+1:end,:) = true;
end

% Seal open ends a bit more: thicken the first/last few axial layers wherever solid starts/ends
axial_skin_vox = max(1, round(0.02 * NzTotal));
axial_skin_vox = min(axial_skin_vox, max(1, floor(NzTotal/2)));
Nt = size(solidMask,1); Nr = size(solidMask,3);
for i = 1:Nt
    for r = 1:Nr
        column = squeeze(solidMask(i,:,r));
        firstSolid = find(column,1,'first');
        if ~isempty(firstSolid) && firstSolid <= axial_skin_vox + 1
            column(1:firstSolid) = true;
        end
        lastSolid = find(column,1,'last');
        if ~isempty(lastSolid) && lastSolid >= (NzTotal - axial_skin_vox)
            column(lastSolid:end) = true;
        end
        solidMask(i,:,r) = column;
    end
end
% Periodic wrap in theta (t = 1 equals t = Nt+1)
solidMask(:,end,:) = solidMask(:,1,:);

solidFractionActual = mean(solidMask(:));
spinodalType = 'anisotropic';
if all(cfg.coneDeg >= 89), spinodalType = 'isotropic'; end

runTimestamp = datetime('now');
runLabelBase = sprintf('stentTiled_%s_Ri%03dum_Ro%03dum_H%03dum_tiles%d', ...
    lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), round(1e6*cfg.H), cfg.tilesAxial);
runDir = unique_run_dir(resultsRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('spinodoid_STENT_TILED_Ri%03dum_Ro%03dum_H%03dum_tiles%d.stl', ...
    round(1e6*cfg.Ri), round(1e6*cfg.Ro), round(1e6*cfg.H), cfg.tilesAxial);
stlPath = unique_path(fullfile(runDir, stlFileName));

% Build coordinate vectors
thetaVec = linspace(0, 2*pi, cfg.Nt+1); thetaVec(end) = [];
if NzTotal > 1, zVec = linspace(0, HTotal, NzTotal); else, zVec = 0; end
if cfg.Nr > 1, rVec = linspace(cfg.Ri, cfg.Ro, cfg.Nr); else, rVec = cfg.Ri; end

% Pad + wrap for isosurface
maskRTZ = permute(solidMask, [3 1 2]); % [Nr Nt Nz] -> [r t z]
[maskRadial, rPad]     = padRadial(maskRTZ, rVec);
[maskAxial,  zPad]     = padAxial(maskRadial, zVec);
[maskWrapped, thetaWrap] = wrapTheta(maskAxial, thetaVec);
[XGrid, YGrid, ZGrid]  = cylindricalGrid(rPad, thetaWrap, zPad);

fv = isosurface(XGrid, YGrid, ZGrid, double(maskWrapped), 0.5);
if isempty(fv.vertices)
    error('tile_spinodoid_stent:emptySurface', 'Stent surface is empty. Check parameters.');
end

meshStats.numFaces = size(fv.faces,1);
meshStats.numVertices = size(fv.vertices,1);

stlwrite_ascii(stlPath, fv.vertices, double(fv.faces));
fprintf('Wrote STL: %s\n', stlPath);

% Log
logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');
elapsedSeconds = toc(tStart);
paramLines = {
    sprintf('  Ri_mm: %.3f', cfg.Ri*1e3)
    sprintf('  Ro_mm: %.3f', cfg.Ro*1e3)
    sprintf('  H_mm (single): %.3f', cfg.H*1e3)
    sprintf('  tiles_axial: %d', cfg.tilesAxial)
    sprintf('  H_total_mm: %.3f', HTotal*1e3)
    sprintf('  Nt/Nz/Nr base: %d / %d / %d', cfg.Nt, cfg.Nz, cfg.Nr)
    sprintf('  lambda_vox: %.3f', cfg.lambda_vox)
    sprintf('  bandwidth: %.3f', cfg.bandwidth)
    sprintf('  nModes: %d', cfg.nModes)
    sprintf('  solid_frac_target: %.3f', cfg.solid_frac)
    sprintf('  coneDeg: [%s]', num2str(cfg.coneDeg))
    sprintf('  sigma_vox: %.3f', cfg.sigma_vox)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    };
geomLines = {
    sprintf('  outer_skin_vox: %d', cfg.outer_skin_vox)
    sprintf('  solid_frac_actual: %.3f', solidFractionActual)
    sprintf('  modes_used: %d', meta.nModes)
    };

logLines = [
    {
    'Spinodoid stent (tiled) run log'
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

fprintf('\n---- TILED STENT SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Ri/Ro/H_total (mm): %.3f / %.3f / %.3f\n', cfg.Ri*1e3, cfg.Ro*1e3, HTotal*1e3);
fprintf('Tiles (axial): %d\n', cfg.tilesAxial);
fprintf('Solid fraction (target/actual): %.3f / %.3f\n', cfg.solid_frac, solidFractionActual);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('STL: %s\n', stlPath);
fprintf('-----------------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
    fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end


% ===================== LOCAL FUNCTIONS =====================

function [solidMask, meta] = generate_stent_mask(cfg)
% Build periodic spinodal field in (theta,z,r) logical grid
Ltheta = 2*pi * ((cfg.Ri + cfg.Ro)/2);
Lz = cfg.H;
Lr = cfg.Ro - cfg.Ri;

[phi, meta] = spinodal_periodic_field_rect([cfg.Nt cfg.Nz cfg.Nr], ...
                                           [Ltheta Lz Lr], ...
                                           cfg.lambda_vox, cfg.bandwidth, ...
                                           cfg.nModes, cfg.coneDeg);

% Optional blur/normalize
if cfg.sigma_vox > 0
    phi = periodic_gaussian_blur(phi, cfg.sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

% Threshold to target solid fraction
t = prctile(phi(:), 100*cfg.solid_frac);
solidMask = phi > t;

% Keep largest component if requested
if cfg.keepLargest && exist('bwconncomp','file') == 2
    try
        CC = bwconncomp(solidMask, 26);
        if CC.NumObjects > 0
            [~,iMax] = max(cellfun(@numel, CC.PixelIdxList));
            tmp = false(size(solidMask));
            tmp(CC.PixelIdxList{iMax}) = true;
            solidMask = tmp;
        end
    catch
        % leave as-is if toolbox not available
    end
end
end

function [maskPad, rPad] = padRadial(mask, rVec)
% mask: [Nr Nt Nz]
dr   = spacing_stent(rVec);
rPad = [max(rVec(1)-dr,0), rVec(:).', rVec(end)+dr];
maskPad = false(size(mask)+[2 0 0]);
maskPad(2:end-1,:,:) = mask;
end

function [maskPad, zPad] = padAxial(mask, zVec)
dz   = spacing_stent(zVec);
zPad = [zVec(1)-dz, zVec(:).', zVec(end)+dz];
maskPad = false(size(mask)+[0 0 2]);
maskPad(:,:,2:end-1) = mask;
end

function [maskWrap, thetaWrap] = wrapTheta(mask, thetaVec)
% Periodic wrap in theta (append first slice to the end)
dtheta   = spacing_stent(thetaVec);
thetaWrap = [thetaVec(:).', thetaVec(end)+dtheta];
maskWrap  = cat(2, mask, mask(:,1,:));
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
