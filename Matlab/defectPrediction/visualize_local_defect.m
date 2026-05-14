% VISUALIZE_LOCAL_DEFECT  Quick sanity check for a single local-defect case.
%
%   Loads sheet.mat and defect_case.json from the given case directory,
%   then plots a top-down view and a mid-z cross-section of the spin layer.
%
% Usage:
%   caseDir = 'defectPrediction/results/lam000/crack_center';
%   run('Matlab/defectPrediction/visualize_local_defect.m')
%
% Or call directly after setting caseDir in the workspace.

if ~exist('caseDir', 'var')
    error('visualize_local_defect:noCaseDir', ...
        'Set caseDir to the case directory before running.');
end

caseDir = local_resolve_case_dir(caseDir);

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

% sheetMask is [row=x_phys, col=y_phys, z]; transpose so image() puts
% physical x on the display x-axis.
topDown = squeeze(any(spinLayer, 3)).';

% Try to load the matching baseline (same angle, same Matlab/...) so we can
% colour the defect (material removed) separately from the spinodoid.
baselineTopDown = local_load_baseline_topdown(caseDir, tbV);
defectMask = false(size(topDown));
if ~isempty(baselineTopDown) && isequal(size(baselineTopDown), size(topDown))
    defectMask = baselineTopDown & ~topDown;
end

% Mid-z cross-section.
midZ     = round(tsV / 2);
midSlice = spinLayer(:,:,max(1, midZ)).';

% Parse metadata title and defect outline.
titleStr = '';
outlines = {};
if isfile(metaPath)
    meta = jsondecode(fileread(metaPath));
    titleStr = sprintf('%s | angle=%d° | type=%s | pos=%s', ...
        meta.case_name, meta.lamellar_angle_deg, meta.defect_type, meta.position_name);
    if isfield(meta, 'local_defects') && ~isempty(meta.local_defects)
        outlines = local_defect_outlines(meta.local_defects);
    end
end

Lx = N * Svox * 1e3;  % mm

fig = figure('Name', titleStr, 'Position', [100 100 900 420]);
subplot(1,2,1);
rgbTop = local_compose_rgb(topDown, defectMask);
image([0 Lx], [0 Lx], rgbTop);
axis equal tight;
set(gca, 'YDir', 'normal');
xlim([0 Lx]); ylim([0 Lx]);
hold on;
local_overlay_outlines(outlines, Lx);
hold off;
xlabel('x (mm)');  ylabel('y (mm)');
title('Spin layer — top-down (blue = material, red = defect)');

subplot(1,2,2);
midDefect = false(size(midSlice));
if ~isempty(baselineTopDown)
    midDefect = defectMask;
end
rgbMid = local_compose_rgb(midSlice, midDefect);
image([0 Lx], [0 Lx], rgbMid);
axis equal tight;
set(gca, 'YDir', 'normal');
xlim([0 Lx]); ylim([0 Lx]);
hold on;
local_overlay_outlines(outlines, Lx);
hold off;
xlabel('x (mm)');  ylabel('y (mm)');
title(sprintf('Spin layer — z-slice %d/%d', midZ, tsV));

if ~isempty(titleStr)
    sgtitle(titleStr, 'Interpreter', 'none');
end

function rgb = local_compose_rgb(solidMask, defectMask)
% Light-blue spinodoid, red defect, white background.
spinColor   = [0.60 0.78 0.95];   % light blue
defectColor = [0.85 0.18 0.18];   % red
bgColor     = [1 1 1];            % white

[H, W] = size(solidMask);
rgb = repmat(reshape(bgColor, 1, 1, 3), H, W);

solidIdx = solidMask & ~defectMask;
defIdx   = defectMask;

for c = 1:3
    chan = rgb(:,:,c);
    chan(solidIdx) = spinColor(c);
    chan(defIdx)   = defectColor(c);
    rgb(:,:,c) = chan;
end
end

function outlines = local_defect_outlines(localDefects)
% Build a cell array of [Nx2] polylines (in mm) describing each defect's full
% intended footprint (independent of where it actually intersects material).
% apply_local_defects writes spinLayer(iy, ix, :), which reflects the intended
% (x,y) across y=x in physical space. Mirror the outline so it aligns with
% the stamped defect in the displayed image.
outlines = {};
if isstruct(localDefects) && isscalar(localDefects)
    specs = localDefects;
elseif isstruct(localDefects)
    specs = localDefects;     % struct array
elseif iscell(localDefects)
    specs = [localDefects{:}];
else
    return;
end
for i = 1:numel(specs)
    d = specs(i);
    if ~isfield(d, 'type') || isempty(d.type), continue; end
    x0_mm = d.position(1) * 1e3;
    y0_mm = d.position(2) * 1e3;
    switch lower(d.type)
        case 'crack'
            L_mm  = d.length_m * 1e3;
            W_mm  = 0;
            if isfield(d, 'width_m'), W_mm = d.width_m * 1e3; end
            W_mm  = max(W_mm, 0.2);   % minimum visible width
            th    = deg2rad(d.theta_deg);
            cx = cos(th); sx = sin(th);
            hx = (L_mm/2) * [cx; sx];
            hw = (W_mm/2) * [-sx; cx];
            corners = [hx + hw, -hx + hw, -hx - hw, hx - hw, hx + hw];
            poly = [x0_mm + corners(1,:); y0_mm + corners(2,:)]';
            outlines{end+1} = poly(:, [2 1]); %#ok<AGROW>
        case 'hole'
            r_mm  = d.radius_m * 1e3;
            count = 1;
            if isfield(d, 'count'), count = d.count; end
            spiral = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; 1 -1; -1 -1];
            th = linspace(0, 2*pi, 96);
            for ci = 1:count
                off = spiral(min(ci, size(spiral,1)), :) * 2.5 * r_mm;
                cxm = x0_mm + off(1);
                cym = y0_mm + off(2);
                poly = [cxm + r_mm * cos(th); cym + r_mm * sin(th)]';
                outlines{end+1} = poly(:, [2 1]); %#ok<AGROW>
            end
    end
end
end

function local_overlay_outlines(outlines, Lx)
% Plot each polyline as a red dotted line, tiled across periodic neighbours so
% outlines that straddle the domain edge wrap correctly.
if isempty(outlines), return; end
offsets = [0 0; Lx 0; -Lx 0; 0 Lx; 0 -Lx; Lx Lx; -Lx Lx; Lx -Lx; -Lx -Lx];
for i = 1:numel(outlines)
    P = outlines{i};
    for k = 1:size(offsets, 1)
        plot(P(:,1) + offsets(k,1), P(:,2) + offsets(k,2), ...
            ':', 'Color', [0.85 0.18 0.18], 'LineWidth', 1.4);
    end
end
end

function topDown = local_load_baseline_topdown(caseDir, tbV)
% Look for sibling 'baseline' directory containing sheet.mat. Returns [] if
% not found or unreadable. Used to derive the defect mask via mask diff.
topDown = [];
parent = fileparts(caseDir);
baseDir = fullfile(parent, 'baseline');
if strcmp(caseDir, baseDir)
    return;   % current case IS baseline; nothing to subtract.
end
matPath = fullfile(baseDir, 'sheet.mat');
if ~isfile(matPath), return; end
try
    Sb = load(matPath, 'sheetMask');
    mb = Sb.sheetMask;
    spin = mb(:,:,tbV+1:end);
    topDown = squeeze(any(spin, 3)).';
catch
    topDown = [];
end
end

function resolved = local_resolve_case_dir(caseDir)
% Resolve caseDir against (1) as-given, (2) repo root, (3) matlab root.
% Lets the script work whether the user is in the repo root or in Matlab/.
if isfolder(caseDir)
    resolved = caseDir;
    return;
end
scriptDir  = fileparts(mfilename('fullpath'));   % .../Matlab/defectPrediction
matlabRoot = fileparts(scriptDir);                % .../Matlab
repoRoot   = fileparts(matlabRoot);               % repo root
candidates = {fullfile(repoRoot, caseDir), fullfile(matlabRoot, caseDir)};
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        resolved = candidates{i};
        return;
    end
end
error('visualize_local_defect:noCaseDir', ...
    'Case directory not found: %s (also tried under repo root and Matlab/)', caseDir);
end
