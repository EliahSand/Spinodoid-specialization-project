function PSS(params)
tStart = tic;
%---PSS - Periodic Spinodal Stent----
%Axially tiles a spinodal stent with angular periodization across theta partitions.

if nargin < 1, params = struct(); end

cfg = struct('Ri',5e-3, ...
    'Ro',6e-3, ...
    'H',20e-3, ...
    'Nt',256, ...
    'Nz',128, ...
    'Nr',48, ...
    'lambda_vox',32, ...
    'bandwidth',0.25, ...
    'nModes',4000, ...
    'solid_frac',0.5, ...
    'coneDeg',[90 90 90], ...
    'rngSeed',1, ...
    'sigma_vox',0.0, ...
    'keepLargest',false, ...
    'outer_skin_vox',0, ...
    'tilesAxial',3, ...
    'theta_partitions',4);
fields = fieldnames(cfg);
for i = 1:numel(fields)
    fn = fields{i};
    if isfield(params, fn) && ~isempty(params.(fn)), cfg.(fn) = params.(fn); end
end

cfg.tilesAxial = max(1, round(cfg.tilesAxial));
cfg.Nt = max(8, round(cfg.Nt));
cfg.Nz = max(4, round(cfg.Nz));
cfg.Nr = max(3, round(cfg.Nr));
cfg.theta_partitions = max(1, round(cfg.theta_partitions));
if cfg.Ro <= cfg.Ri, error('PSS:invalidRadius', 'Ro must exceed Ri'); end

scriptDir = fileparts(mfilename('fullpath'));
helpersDir = fullfile(scriptDir,'helpers');
if ~isfolder(helpersDir), error('Helpers directory missing: %s', helpersDir); end
addpath(helpersDir);

%% Generate base stent mask with angular periodization
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

phi = enforce_theta_partitions(phi, cfg.theta_partitions);
cfg.Nt = size(phi,1); % update if partitions changed size

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
        % keep full mask if toolbox unavailable
    end
end

solidMask(:,end,:) = solidMask(:,1,:); % theta periodicity

solidFractionSegment = mean(solidMask(:));

%% Tile along axial direction
solidMaskTiled = repmat(solidMask, [1 cfg.tilesAxial 1]);
NzTotal = size(solidMaskTiled,2);
HTotal = cfg.H * cfg.tilesAxial;

if cfg.outer_skin_vox > 0
    s = min(cfg.outer_skin_vox, floor(NzTotal/8));
    solidMaskTiled(:,1:s,:) = true;
    solidMaskTiled(:,end-s+1:end,:) = true;
end

axial_skin_vox_total = max(1, round(0.02 * NzTotal));
axial_skin_vox_total = min(axial_skin_vox_total, max(1, floor(NzTotal/2)));
Nt = size(solidMaskTiled,1);
Nr = size(solidMaskTiled,3);
for i = 1:Nt
    for r = 1:Nr
        column = squeeze(solidMaskTiled(i,:,r));
        firstSolid = find(column, 1, 'first');
        if ~isempty(firstSolid) && firstSolid <= axial_skin_vox_total + 1
            column(1:firstSolid) = true;
        end
        lastSolid = find(column, 1, 'last');
        if ~isempty(lastSolid) && lastSolid >= (NzTotal - axial_skin_vox_total)
            column(lastSolid:end) = true;
        end
        solidMaskTiled(i,:,r) = column;
    end
end

solidMaskTiled(:,end,:) = solidMaskTiled(:,1,:);
solidFractionTiled = mean(solidMaskTiled(:));

%% Build run folder and export STL
spinodalType = 'anisotropic';
if all(cfg.coneDeg >= 89)
    spinodalType = 'isotropic';
elseif cfg.coneDeg(1) == 30 && all(cfg.coneDeg([2 3]) == 0)
    spinodalType = 'lamellar';
elseif cfg.coneDeg(2) == 30 && all(cfg.coneDeg([1 3]) == 0)
    spinodalType = 'columnar';
elseif all(cfg.coneDeg == 15)
    spinodalType = 'cubic';
end

isPeriodic = (cfg.theta_partitions > 1) || (cfg.tilesAxial > 1);
if isPeriodic
    modeLabel = 'periodic';
else
    modeLabel = 'single';
end
resultsRoot = fullfile(scriptDir,'results','stents',modeLabel);
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

runTimestamp = datetime('now');
if isPeriodic
    runLabelBase = sprintf('stent_%s_%s_Ri%03dum_Ro%03dum_H%03dum_tp%d_tiles%d', ...
        modeLabel, lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), ...
        round(1e6*cfg.H), cfg.theta_partitions, cfg.tilesAxial);
else
    runLabelBase = sprintf('stent_%s_%s_Ri%03dum_Ro%03dum_H%03dum', ...
        modeLabel, lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), ...
        round(1e6*cfg.H));
end
typeRoot = fullfile(resultsRoot, lower(spinodalType));
if ~exist(typeRoot,'dir'), mkdir(typeRoot); end
runDir = unique_run_dir(typeRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlPath = unique_path(fullfile(runDir, 'stent.stl'));
[~, stlFileName, ext] = fileparts(stlPath);
stlFileName = [stlFileName ext];

thetaVec = linspace(0, 2*pi, Nt+1); thetaVec(end) = [];
if NzTotal > 1, zVec = linspace(0, HTotal, NzTotal); else, zVec = 0; end
if cfg.Nr > 1, rVec = linspace(cfg.Ri, cfg.Ro, cfg.Nr); else, rVec = cfg.Ri; end

maskRTZ = permute(solidMaskTiled, [3 1 2]); % [Nr Nt Nz] -> [r t z]
[maskRadial, rPad]     = padRadial_stent(maskRTZ, rVec);
[maskAxial,  zPad]     = padAxial_stent(maskRadial, zVec);
[maskWrapped, thetaWrap] = wrapTheta_stent(maskAxial, thetaVec);
[XGrid, YGrid, ZGrid]  = cylindricalGrid_stent(rPad, thetaWrap, zPad);

fv = isosurface(XGrid, YGrid, ZGrid, double(maskWrapped), 0.5);
if isempty(fv.vertices)
    error('PSS:emptySurface', 'Stent surface is empty. Check parameters.');
end

fv.faces = double(fv.faces);
meshStats.numFaces = size(fv.faces,1);
meshStats.numVertices = size(fv.vertices,1);

stlwrite_ascii(stlPath, fv.vertices, double(fv.faces));
fprintf('Wrote STL: %s\n', stlPath);

[matDir, matBase] = fileparts(stlPath);
matPath = fullfile(matDir, sprintf('%s.mat', matBase));
maskTheta = thetaVec;
maskZ = zVec;
maskR = rVec;
save(matPath, 'solidMaskTiled', 'maskTheta', 'maskZ', 'maskR', '-v7.3');
fprintf('Saved voxel mask to: %s\n', matPath);

thetaSpacing = 2*pi * ((cfg.Ri + cfg.Ro)/2) / Nt;
axialSpacing = HTotal / NzTotal;
radialSpacing = (cfg.Ro - cfg.Ri) / max(1, cfg.Nr - 1);
manifest = struct();
manifest.mask = [matBase '.mat'];
manifest.var = 'solidMaskTiled';
manifest.spacing = [thetaSpacing axialSpacing radialSpacing];
manifest.origin = [0 0 0];
manifest.material = 'SPINODAL';
manifest.notes = 'Axes correspond to [theta, z, r] in cylindrical coordinates.';
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

%% Logging
logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');
elapsedSeconds = toc(tStart);

paramLines = {
    sprintf('  Ri_mm: %.3f', cfg.Ri*1e3)
    sprintf('  Ro_mm: %.3f', cfg.Ro*1e3)
    sprintf('  H_mm (segment): %.3f', cfg.H*1e3)
    sprintf('  tiles_axial: %d', cfg.tilesAxial)
    sprintf('  H_total_mm: %.3f', HTotal*1e3)
    sprintf('  Nt/Nz/Nr (segment): %d / %d / %d', Nt, cfg.Nz, cfg.Nr)
    sprintf('  lambda_vox: %.3f', cfg.lambda_vox)
    sprintf('  bandwidth: %.3f', cfg.bandwidth)
    sprintf('  nModes: %d', cfg.nModes)
    sprintf('  solid_frac_target: %.3f', cfg.solid_frac)
    sprintf('  coneDeg: [%s]', num2str(cfg.coneDeg))
    sprintf('  sigma_vox: %.3f', cfg.sigma_vox)
    sprintf('  theta_partitions: %d', cfg.theta_partitions)
    sprintf('  mode: %s', modeLabel)
    sprintf('  keep_largest: %d', cfg.keepLargest)
    sprintf('  outer_skin_vox: %d', cfg.outer_skin_vox)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    sprintf('  modes_used: %d', meta.nModes)
    };

geomLines = {
    sprintf('  solid_frac_segment: %.3f', solidFractionSegment)
    sprintf('  solid_frac_tiled: %.3f', solidFractionTiled)
    };

logLines = [
    {
    'Spinodoid stent (angular periodized, tiled) run log'
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

fprintf('\n---- PERIODIZED TILED STENT SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Ri/Ro/H (mm): %.3f / %.3f / %.3f\n', cfg.Ri*1e3, cfg.Ro*1e3, cfg.H*1e3);
fprintf('Theta partitions: %d (%.1f deg sectors)\n', cfg.theta_partitions, 360/cfg.theta_partitions);
fprintf('Tiles (axial): %d -> total H %.3f mm\n', cfg.tilesAxial, HTotal*1e3);
fprintf('Grid Nt/Nz/Nr (segment): %d / %d / %d\n', Nt, cfg.Nz, cfg.Nr);
fprintf('Solid fraction (segment/tiled): %.3f / %.3f\n', solidFractionSegment, solidFractionTiled);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('STL path: %s\n', stlPath);
fprintf('----------------------------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
    fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end

function phiRep = enforce_theta_partitions(phi, thetaPartitions)
phiRep = phi;
Nt = size(phi,1);
if thetaPartitions <= 1, return; end
if mod(Nt, thetaPartitions) ~= 0
    error('PSS:thetaDiv', ...
        'Nt (%d) must be divisible by theta_partitions (%d).', Nt, thetaPartitions);
end
NtSector = Nt / thetaPartitions;
phiReshaped = reshape(phi, NtSector, thetaPartitions, size(phi,2), size(phi,3));
sectorBlend = mean(phiReshaped, 2);
sectorBlend = reshape(sectorBlend, NtSector, 1, size(phi,2), size(phi,3));
phiRep = repmat(sectorBlend, 1, thetaPartitions, 1, 1);
phiRep = reshape(phiRep, Nt, size(phi,2), size(phi,3));
end

function [maskPad, rPad] = padRadial_stent(mask, rVec)
dr = spacing_stent(rVec);
rPad = [max(rVec(1) - dr, 0), rVec(:).', rVec(end) + dr];
maskPad = false(size(mask) + [2 0 0]);
maskPad(2:end-1,:,:) = mask;
end

function [maskPad, zPad] = padAxial_stent(mask, zVec)
dz = spacing_stent(zVec);
zPad = [zVec(1) - dz, zVec(:).', zVec(end) + dz];
maskPad = false(size(mask) + [0 0 2]);
maskPad(:,:,2:end-1) = mask;
end

function [maskWrap, thetaWrap] = wrapTheta_stent(mask, thetaVec)
dtheta = spacing_stent(thetaVec);
thetaWrap = [thetaVec(:).', thetaVec(end) + dtheta];
maskWrap = cat(2, mask, mask(:,1,:));
end

function [X, Y, Z] = cylindricalGrid_stent(rVec, thetaVec, zVec)
[R, Theta, Z] = ndgrid(rVec, thetaVec, zVec);
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
