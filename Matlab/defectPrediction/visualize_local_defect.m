% VISUALIZE_LOCAL_DEFECT  Quick sanity check for a single local-defect case.
%
%   Loads sheet.mat and defect_case.json from the given case directory,
%   then plots a top-down view and a mid-z cross-section of the spin layer.
%
% Usage:
%   caseDir = 'Matlab/defectPrediction/results/lam000/crack_center';
%   run('Matlab/defectPrediction/visualize_local_defect.m')
%
% Or call directly after setting caseDir in the workspace.

if ~exist('caseDir', 'var')
    error('visualize_local_defect:noCaseDir', ...
        'Set caseDir to the case directory before running.');
end

matPath  = fullfile(caseDir, 'sheet.mat');
metaPath = fullfile(caseDir, 'defect_case.json');

if ~isfile(matPath)
    error('visualize_local_defect:noMat', 'sheet.mat not found in %s', caseDir);
end

S         = load(matPath, 'sheetMask', 'voxelSizeXY', 'zVoxelThickness');
sheetMask = S.sheetMask;
Svox      = S.voxelSizeXY;

% Reconstruct layer boundaries from zVoxelThickness.
dz       = S.zVoxelThickness;
Nz       = numel(dz);
% Assume first layers are base (thicker dz) and last are spin (thinner dz).
% Use PSSCone convention: t_base=2e-3, t_spin=1e-3; split at tbV.
baseParams_L = 40e-3;
N  = size(sheetMask, 1);
t_base = 2e-3;  t_spin = 1e-3;
tbV = max(3, round(t_base / (baseParams_L / N)));
tsV = Nz - tbV;

spinLayer  = sheetMask(:,:,tbV+1:end);

% Top-down projection (solid if any z in spin layer is solid).
topDown = squeeze(any(spinLayer, 3));

% Mid-z cross-section.
midZ     = round(tsV / 2);
midSlice = spinLayer(:,:,max(1, midZ));

% Parse metadata title.
titleStr = '';
if isfile(metaPath)
    meta = jsondecode(fileread(metaPath));
    titleStr = sprintf('%s | angle=%d° | type=%s | pos=%s', ...
        meta.case_name, meta.lamellar_angle_deg, meta.defect_type, meta.position_name);
end

Lx = N * Svox * 1e3;  % mm

fig = figure('Name', titleStr, 'Position', [100 100 900 420]);
subplot(1,2,1);
imagesc([0 Lx], [0 Lx], topDown);
colormap(gca, gray);
axis equal tight;
xlabel('x (mm)');  ylabel('y (mm)');
title('Spin layer — top-down projection');

subplot(1,2,2);
imagesc([0 Lx], [0 Lx], midSlice);
colormap(gca, gray);
axis equal tight;
xlabel('x (mm)');  ylabel('y (mm)');
title(sprintf('Spin layer — z-slice %d/%d', midZ, tsV));

if ~isempty(titleStr)
    sgtitle(titleStr, 'Interpreter', 'none');
end
