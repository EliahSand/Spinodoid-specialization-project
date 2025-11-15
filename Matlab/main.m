function main()
tStart = tic;
% ------------------------- Design knobs ---------------------------------------
N            = 128;        % grid size (NxNxN). Use powers of two for speed
L            = 0.5;        % physical box length (mm)
lambda_vox   = 32;         % target feature wavelength in voxels (~rib/ligament spacing)
bandwidth    = 0.50;       % relative shell thickness around target |k| (0.1–0.3)
nModes       = 4000;       % number of Fourier modes to sample (1k–10k typical)
solid_frac   = 0.5;        % volume fraction of SOLID after threshold (0..1)
coneDeg      = [90 90 90]; % cone half-angles about x,y,z (90= isotropic). e.g. [90 90 90]
rngSeed      = 1;          % reproducible
rng(rngSeed);
scriptDir    = fileparts(mfilename('fullpath'));
helpersDir   = fullfile(scriptDir, 'helpers');
if ~isfolder(helpersDir)
    error('Helpers directory missing: %s', helpersDir);
end
addpath(helpersDir);
resultsRoot = fullfile(scriptDir, 'results', 'cells');
if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end

% ---- Generate periodic field -------------------------------------------
[phi, meta]  = spinodal_periodic_field(N, L, lambda_vox, bandwidth, nModes, coneDeg);

% ---- Threshold by percentile to hit the exact solid fraction ------------
t            = prctile(phi(:), 100*solid_frac);   %note: solid NOT(1-solid_frac)
solidMask    = phi > t;

% ---- Quick sanity prints ------------------------------------------------
fprintf('Field mean=%.3g, std=%.3g, threshold t=%.3g, solids=%.2f%% voxels\n', ...
         mean(phi(:)), std(phi(:)), t, 100*mean(solidMask(:)));

% ---- Make exterior CLOSED (solid skin around the box) --------------------
skin_thickness_vox = max(1, round(0.02*N));   % ~2% of box size; adjust as needed

if skin_thickness_vox > 0
    s = skin_thickness_vox;
    solidMask(1:s,:,:)                 = true;  % -X face
    solidMask(end-s+1:end,:,:)         = true;  % +X face
    solidMask(:,1:s,:)                 = true;  % -Y face
    solidMask(:,end-s+1:end,:)         = true;  % +Y face
    solidMask(:,:,1:s)                 = true;  % -Z face
    solidMask(:,:,end-s+1:end)         = true;  % +Z face
end

% (Optional) Keep only the largest solid component, remove isolated specks
try
    if exist('bwconncomp','file') == 2
        CC = bwconncomp(solidMask, 26);
        [~,iMax] = max(cellfun(@numel, CC.PixelIdxList));
        tmp = false(size(solidMask));
        tmp(CC.PixelIdxList{iMax}) = true;
        solidMask = tmp;
    else
        % Image Processing Toolbox not available — skip connectivity pruning
    end
catch
    % If any error occurs (e.g., toolbox missing), leave solidMask as-is
end

solidFractionActual = mean(solidMask(:));

% ---- Export watertight SOLID surface + log ------------------------------

spinodalType = 'anisotropic';
if all(coneDeg >= 89)
    spinodalType = 'isotropic';
elseif coneDeg(1) == 30 && all(coneDeg([2 3]) == 0)
    spinodalType = 'lamellar';
elseif coneDeg(2) == 30 && all(coneDeg([1 3]) == 0)
    spinodalType = 'columnar';

end

runTimestamp = datetime('now');
runLabelBase = sprintf('%s_N%d_sf%02d', lower(spinodalType), N, round(100*solid_frac));
typeRoot = fullfile(resultsRoot, lower(spinodalType));
if ~exist(typeRoot, 'dir'), mkdir(typeRoot); end
runDir = fullfile(typeRoot, runLabelBase);
suffix = 2;
while exist(runDir, 'dir')
    runDir = fullfile(typeRoot, sprintf('%s_run%02d', runLabelBase, suffix));
    suffix = suffix + 1;
end
if ~exist(runDir, 'dir'), mkdir(runDir); end
[~, runFolderName] = fileparts(runDir);

stlFileName = sprintf('spinodoid_N%d_sf%02d_vox.stl', N, round(100*solid_frac));
stlPath = fullfile(runDir, stlFileName);
if isfile(stlPath)
    [folder, base, ext] = fileparts(stlPath);
    v = 1;
    while isfile(fullfile(folder, sprintf('%s_v%d%s', base, v, ext)))
        v = v + 1;
    end
    stlPath = fullfile(folder, sprintf('%s_v%d%s', base, v, ext));
    [~, stlFileName, ext] = fileparts(stlPath);
    stlFileName = [stlFileName ext];
end

meshStats = binary_to_watertight_stl(solidMask, L, stlPath);
fprintf('Wrote file to: %s\n', stlPath);

logPath = fullfile(runDir, 'run_log.txt');
timestampStr = datestr(runTimestamp, 'dd-mmm-yyyy HH:MM:SS');
elapsedSeconds = toc(tStart);

paramLines = {
    sprintf('  grid: %d', N)
    sprintf('  relative_density: %.3f', solid_frac)
    sprintf('  box_length: %.3f', L)
    sprintf('  lambda_vox: %.3f', lambda_vox)
    sprintf('  bandwidth: %.3f', bandwidth)
    sprintf('  nModes: %d', nModes)
    sprintf('  coneDeg: [%s]', num2str(coneDeg))
    sprintf('  skin_thickness_vox: %d', skin_thickness_vox)
    sprintf('  rng_seed: %d', rngSeed)
    };

fieldMetaLines = {
    sprintf('  N: %d', meta.N)
    sprintf('  L: %.3f', meta.L)
    sprintf('  lambda_vox: %.3f', meta.lambda_vox)
    sprintf('  k_target_idx: %.3f', meta.k_target_idx)
    sprintf('  bandwidth: %.3f', meta.bandwidth)
    sprintf('  nModes: %d', meta.nModes)
    };

logLines = [
    {
    'Spinodal run log'
    sprintf('Generated: %s', timestampStr)
    sprintf('Run folder: %s', runFolderName)
    sprintf('Type: %s', spinodalType)
    'Mode: periodic spinodal'
    sprintf('STL: %s', stlFileName)
    sprintf('Faces: %d', meshStats.numFaces)
    sprintf('Vertices: %d', meshStats.numVertices)
    sprintf('Volume fraction (target/actual): %.3f / %.3f', solid_frac, solidFractionActual)
    ''
    'Parameters:'
    } ;
    paramLines ;
    {
    ''
    'Field metadata:'
    } ;
    fieldMetaLines
    {
    ''
    'Runtime:'
    sprintf('  elapsed_seconds: %.2f', elapsedSeconds)
    sprintf('  elapsed_minutes: %.2f', elapsedSeconds/60)
    }
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
fprintf('Type/Mode: %s / periodic spinodal\n', spinodalType);
fprintf('Grid: %d^3, lambda_vox: %.3f, bandwidth: %.3f\n', N, lambda_vox, bandwidth);
fprintf('Volume fraction (target/actual): %.3f / %.3f\n', solid_frac, solidFractionActual);
fprintf('STL: %s\n', stlPath);
fprintf('Faces: %d, Vertices: %d\n', meshStats.numFaces, meshStats.numVertices);
fprintf('Cone degrees: [%s]\n', num2str(coneDeg));
fprintf('----------------------\n\n');
if elapsedSeconds < 60
    fprintf('Run completed in %.2f seconds.\n', elapsedSeconds);
else
    fprintf('Run completed in %.2f minutes.\n', elapsedSeconds/60);
end
end
