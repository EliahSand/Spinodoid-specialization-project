function PSSL(params)
tStart = tic;
%----PSSL - Periodic Spinodal Stent Layered----
%Axially tiled relief stent with angularly periodized spinodal layers.

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
    'outer_skin_vox', 0, ...
    'slice_count', 10, ...
    't_spin', 0.5e-3, ...
    't_base', 0.5e-3, ...
    'theta_partitions', 2, ...
    'tilesAxial', 2);

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
cfg.tilesAxial = max(1, round(cfg.tilesAxial));
cfg.slice_count = max(1, round(cfg.slice_count));
cfg.outer_skin_vox = max(0, round(cfg.outer_skin_vox));
cfg.theta_partitions = max(1, round(cfg.theta_partitions));

scriptDir  = fileparts(mfilename('fullpath'));
helpersDir = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);

%% Geometry generation (mirrors tile_spinodoid_relief_stent_periodized)
[geomState, meta] = generate_relief_stent_geometry(cfg);
stentMask                 = geomState.stentMask;
Nt                        = geomState.Nt;
NzTotal                   = geomState.NzTotal;
HTotal                    = geomState.HTotal;
shellThickness            = geomState.shellThickness;
nrBase                    = geomState.nrBase;
nrSpin                    = geomState.nrSpin;
actualBaseThickness       = geomState.actualBaseThickness;
actualSpinThickness       = geomState.actualSpinThickness;
solidFractionBase         = geomState.solidFractionBase;
solidFractionSpinLayer    = geomState.solidFractionSpinLayer;
solidFractionSegment      = geomState.solidFractionSegment;
solidFractionTotal        = geomState.solidFractionTotal;
solidFractionPattern      = geomState.solidFractionPattern;

%% Export tiled geometry
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
resultsRoot = fullfile(scriptDir, 'results', 'stent_relief', modeLabel);
if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

runTimestamp = datetime('now');
if isPeriodic
    runLabelBase = sprintf('stentRelief_%s_%s_Ri%03dum_Ro%03dum_H%03dum_tp%d_tiles%d', ...
        modeLabel, lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), ...
        round(1e6*cfg.H), cfg.theta_partitions, cfg.tilesAxial);
else
    runLabelBase = sprintf('stentRelief_%s_%s_Ri%03dum_Ro%03dum_H%03dum', ...
        modeLabel, lower(spinodalType), round(1e6*cfg.Ri), round(1e6*cfg.Ro), ...
        round(1e6*cfg.H));
end
typeRoot = fullfile(resultsRoot, lower(spinodalType));
if ~exist(typeRoot,'dir'), mkdir(typeRoot); end
runDir = unique_run_dir(typeRoot, runLabelBase);
[~, runFolderName] = fileparts(runDir);

stlPath = unique_path(fullfile(runDir, 'stent_relief.stl'));
[~, stlFileName, ext] = fileparts(stlPath);
stlFileName = [stlFileName ext];

thetaVec = linspace(0, 2*pi, Nt+1); thetaVec(end) = [];
if NzTotal > 1
    zVec = linspace(0, HTotal, NzTotal);
else
    zVec = 0;
end
if cfg.Nr > 1
    rVec = linspace(cfg.Ri, cfg.Ro, cfg.Nr);
else
    rVec = cfg.Ri;
end

maskRTZ = permute(stentMask, [3 1 2]);
[maskRadial, rPad] = padRadial_stent(maskRTZ, rVec);
[maskAxial, zPad] = padAxial_stent(maskRadial, zVec);
[maskWrapped, thetaWrap] = wrapTheta_stent(maskAxial, thetaVec);
[XGrid, YGrid, ZGrid] = cylindricalGrid_stent(rPad, thetaWrap, zPad);

fv = isosurface(XGrid, YGrid, ZGrid, double(maskWrapped), 0.5);
if isempty(fv.vertices)
    error('PSSL:emptySurface', ...
        'Relief stent surface is empty. Adjust parameters.');
end

fv.faces = double(fv.faces);
normals = faceNormals(fv.vertices, fv.faces);
centroids = (fv.vertices(fv.faces(:,1),:) + fv.vertices(fv.faces(:,2),:) + fv.vertices(fv.faces(:,3),:))/3;
radial = sum(centroids(:,1:2).^2, 2);
flipIdx = radial < ((cfg.Ri + cfg.Ro)/2)^2 & dot(normals, centroids, 2) > 0;
fv.faces(flipIdx,:) = fv.faces(flipIdx,[1 3 2]);

meshStats.numFaces = size(fv.faces,1);
meshStats.numVertices = size(fv.vertices,1);

stlwrite_ascii(stlPath, fv.vertices, fv.faces);
fprintf('Wrote STL: %s\n', stlPath);

[matDir, matBase] = fileparts(stlPath);
matPath = fullfile(matDir, sprintf('%s.mat', matBase));
maskTheta = thetaVec;
maskZ = zVec;
maskR = rVec;
save(matPath, 'stentMask', 'maskTheta', 'maskZ', 'maskR', '-v7.3');
fprintf('Saved voxel mask to: %s\n', matPath);

thetaSpacing = 2*pi * ((cfg.Ri + cfg.Ro)/2) / Nt;
axialSpacing = HTotal / NzTotal;
radialSpacing = (cfg.Ro - cfg.Ri) / max(1, cfg.Nr - 1);
manifest = struct();
manifest.mask = [matBase '.mat'];
manifest.var = 'stentMask';
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
    sprintf('  lambda_vox: %.3f', cfg.lambda_vox)
    sprintf('  bandwidth: %.3f', cfg.bandwidth)
    sprintf('  nModes: %d', cfg.nModes)
    sprintf('  solid_frac_target: %.3f', cfg.solid_frac)
    sprintf('  coneDeg: [%s]', num2str(cfg.coneDeg))
    sprintf('  sigma_vox: %.3f', cfg.sigma_vox)
    sprintf('  slice_count: %d', cfg.slice_count)
    sprintf('  mode: %s', modeLabel)
    sprintf('  theta_partitions: %d', cfg.theta_partitions)
    sprintf('  rng_seed: %d', cfg.rngSeed)
    sprintf('  modes_used: %d', meta.nModes)
    };

stentLines = {
    sprintf('  Ri_mm: %.3f', cfg.Ri*1e3)
    sprintf('  Ro_mm: %.3f', cfg.Ro*1e3)
    sprintf('  H_mm (segment): %.3f', cfg.H*1e3)
    sprintf('  tiles_axial: %d', cfg.tilesAxial)
    sprintf('  H_total_mm: %.3f', HTotal*1e3)
    sprintf('  shell_thickness_mm: %.3f', shellThickness*1e3)
    sprintf('  Nt/Nz/Nr (segment): %d / %d / %d', Nt, cfg.Nz, cfg.Nr)
    sprintf('  outer_skin_vox: %d', cfg.outer_skin_vox)
    sprintf('  solid_frac segment/tiled: %.3f / %.3f', solidFractionSegment, solidFractionTotal)
    };

reliefLines = {
    sprintf('  t_base_mm (target/actual): %.4f / %.4f', cfg.t_base*1e3, actualBaseThickness*1e3)
    sprintf('  t_spin_mm (target/actual): %.4f / %.4f', cfg.t_spin*1e3, actualSpinThickness*1e3)
    sprintf('  voxels base/spin: %d / %d', nrBase, nrSpin)
    sprintf('  solid_frac base/spin layer: %.3f / %.3f', solidFractionBase, solidFractionSpinLayer)
    sprintf('  pattern_coverage_2d: %.3f', solidFractionPattern)
    };

logLines = [
    {
    'Spinodoid relief stent (angular periodized, tiled) run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Folder: %s', runFolderName)
    sprintf('Type: %s', spinodalType)
    sprintf('STL: %s', stlFileName)
    sprintf('Faces: %d', meshStats.numFaces)
    sprintf('Vertices: %d', meshStats.numVertices)
    ''
    'Base cell parameters:'
    } ;
    paramLines ;
    {
    ''
    'Stent geometry:'
    } ;
    stentLines ;
    {
    ''
    'Layering metadata:'
    } ;
    reliefLines ;
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

fprintf('\n---- TILED ANGULAR RELIEF STENT SUMMARY ----\n');
fprintf('Folder: %s\n', runFolderName);
fprintf('Type: %s\n', spinodalType);
fprintf('Ri/Ro (mm): %.3f / %.3f\n', cfg.Ri*1e3, cfg.Ro*1e3);
fprintf('Segment H / total H (mm): %.3f / %.3f (tiles=%d)\n', cfg.H*1e3, HTotal*1e3, cfg.tilesAxial);
fprintf('Theta partitions: %d (%.1f deg sectors)\n', cfg.theta_partitions, 360/cfg.theta_partitions);
fprintf('Grid Nt/Nz/Nr (segment): %d / %d / %d\n', Nt, cfg.Nz, cfg.Nr);
fprintf('Base/spin thickness (mm): %.3f / %.3f\n', actualBaseThickness*1e3, actualSpinThickness*1e3);
fprintf('Spin pattern coverage (2D / segment / tiled): %.3f / %.3f / %.3f\n', ...
    solidFractionPattern, solidFractionSegment, solidFractionTotal);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('STL path: %s\n', stlPath);
fprintf('--------------------------------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end

function [state, meta] = generate_relief_stent_geometry(cfg)
shellThickness = cfg.Ro - cfg.Ri;
if shellThickness <= 0
    error('PSSL:invalidRadius', ...
        'Outer radius must exceed inner radius.');
end
if cfg.t_base <= 0 || cfg.t_spin <= 0
    error('PSSL:invalidThickness', ...
        'Both t_base and t_spin must be positive.');
end
layerSum = cfg.t_base + cfg.t_spin;
tolThickness = max(1e-9, 1e-3 * shellThickness);
if abs(layerSum - shellThickness) > tolThickness
    error('PSSL:thicknessMismatch', ...
        't_base + t_spin (%.6f) must equal shell thickness Ro-Ri (%.6f).', ...
        layerSum, shellThickness);
end

targetBaseVox = cfg.Nr * (cfg.t_base / shellThickness);
nrBase = max(1, min(cfg.Nr-1, round(targetBaseVox)));
nrSpin = cfg.Nr - nrBase;
if nrSpin < 1
    error('PSSL:spinVoxels', ...
        'Spin layer collapsed to zero voxels; increase Nr or reduce t_base.');
end

actualBaseThickness = shellThickness * (nrBase / cfg.Nr);
actualSpinThickness = shellThickness - actualBaseThickness;

rng(cfg.rngSeed);
Ltheta = 2*pi * ((cfg.Ri + cfg.Ro)/2);
Lz = cfg.H;
Lr = shellThickness;
[phi, meta] = spinodal_periodic_field_rect([cfg.Nt cfg.Nz cfg.Nr], ...
    [Ltheta Lz Lr], cfg.lambda_vox, cfg.bandwidth, cfg.nModes, cfg.coneDeg);

if cfg.sigma_vox > 0
    phi = periodic_gaussian_blur(phi, cfg.sigma_vox);
    phi = phi - mean(phi(:));
    phi = phi ./ (std(phi(:)) + eps);
end

phi = enforce_theta_partitions(phi, cfg.theta_partitions);
Nt = size(phi,1);

slab = min(max(1, cfg.slice_count), nrSpin);
phiOuter = phi(:,:,end-slab+1:end);
phi2 = mean(phiOuter, 3);
phi2 = phi2 - mean(phi2(:));
phi2 = phi2 ./ (std(phi2(:)) + eps);

threshold2 = prctile(phi2(:), 100 * cfg.solid_frac);
mask2 = phi2 > threshold2;
solidFractionPattern = mean(mask2(:));

baseLayer = true(Nt, cfg.Nz, nrBase);
spinLayer = repmat(mask2, 1, 1, nrSpin);
segmentMask = cat(3, baseLayer, spinLayer);

solidFractionBase = mean(baseLayer(:));
solidFractionSpinLayer = mean(spinLayer(:));
solidFractionSegment = mean(segmentMask(:));

stentMask = repmat(segmentMask, [1 cfg.tilesAxial 1]);
NzTotal = size(stentMask, 2);
HTotal = cfg.H * cfg.tilesAxial;

if cfg.outer_skin_vox > 0
    s = min(cfg.outer_skin_vox, floor(NzTotal/8));
    stentMask(:,1:s,:) = true;
    stentMask(:,end-s+1:end,:) = true;
end

axial_skin_vox = max(1, round(0.02 * NzTotal));
axial_skin_vox = min(axial_skin_vox, max(1, floor(NzTotal/2)));
Nr = size(stentMask,3);
for i = 1:Nt
    for r = 1:Nr
        column = squeeze(stentMask(i,:,r));
        firstSolid = find(column, 1, 'first');
        if ~isempty(firstSolid) && firstSolid <= axial_skin_vox + 1
            column(1:firstSolid) = true;
        end
        lastSolid = find(column, 1, 'last');
        if ~isempty(lastSolid) && lastSolid >= (NzTotal - axial_skin_vox)
            column(lastSolid:end) = true;
        end
        stentMask(i,:,r) = column;
    end
end

stentMask(:,end,:) = stentMask(:,1,:);
solidFractionTotal = mean(stentMask(:));

state = struct( ...
    'stentMask', stentMask, ...
    'Nt', Nt, ...
    'NzTotal', NzTotal, ...
    'HTotal', HTotal, ...
    'shellThickness', shellThickness, ...
    'nrBase', nrBase, ...
    'nrSpin', nrSpin, ...
    'actualBaseThickness', actualBaseThickness, ...
    'actualSpinThickness', actualSpinThickness, ...
    'solidFractionBase', solidFractionBase, ...
    'solidFractionSpinLayer', solidFractionSpinLayer, ...
    'solidFractionSegment', solidFractionSegment, ...
    'solidFractionTotal', solidFractionTotal, ...
    'solidFractionPattern', solidFractionPattern);
end

function phiRep = enforce_theta_partitions(phi, thetaPartitions)
phiRep = phi;
Nt = size(phi,1);
if thetaPartitions <= 1
    return;
end
if mod(Nt, thetaPartitions) ~= 0
    error('PSSL:thetaDiv', ...
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
    d = max(vec, 1);
end
end
