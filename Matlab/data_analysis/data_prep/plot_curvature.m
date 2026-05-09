function plot_curvature(solid, shell, trLabel, thetaDeg, outDir, opts)
%PLOT_CURVATURE Plot mid-plane YZ curvature using baseline nodes.
%   Baseline = original (Y,Z). Solid/Shell are plotted on top using
%   deformed coordinates (Y+U2, Z+U3) by default.

arguments
    solid table
    shell table
    trLabel (1, :) char
    thetaDeg (1, 1) double
    outDir (1, :) char
    opts.UseU2 (1, 1) logical = true
    opts.UseU3 (1, 1) logical = true
    opts.CurvatureCheck struct = struct()
end

ensure_dir(outDir);

solid = solid(:, :);
shell = shell(:, :);

% Baseline coordinates (reference)
Y0 = double(solid.Y);
Z0 = double(solid.Z);

% Deformed coordinates
Ysol = Y0;
Ysh = Y0;
if opts.UseU2
    Ysol = Y0 + double(solid.U2);
    Ysh = Y0 + double(shell.U2);
end

Zsol = Z0;
Zsh = Z0;
if opts.UseU3
    Zsol = Z0 + double(solid.U3);
    Zsh = Z0 + double(shell.U3);
end

% Sort by Y for a smooth profile
[Y0s, sortIdx] = sort(Y0);
Z0s = Z0(sortIdx);
Ysol = Ysol(sortIdx);
Zsol = Zsol(sortIdx);
Ysh = Ysh(sortIdx);
Zsh = Zsh(sortIdx);

ratioLabel = format_ratio_label(trLabel);
thetaLabel = format_theta_label(thetaDeg);

fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
plot(ax, Y0s, Z0s, '-k', 'LineWidth', 1.5, 'DisplayName', 'Baseline');
hold(ax, 'on');
plot(ax, Ysol, Zsol, '-r', 'LineWidth', 1.5, 'DisplayName', 'Solid');
plot(ax, Ysh, Zsh, '-b', 'LineWidth', 1.5, 'DisplayName', 'Shell');
if isfield(opts.CurvatureCheck, 'theory')
    [yTheory, zTheory] = build_theory_arc(Y0s, Z0s, Ysol, Ysh, Zsol, opts.CurvatureCheck.theory);
    if ~isempty(yTheory)
        plot(ax, yTheory, zTheory, '--g', 'LineWidth', 1.5, 'DisplayName', 'Theory');
    end
end
grid(ax, 'on');
axis(ax, 'equal');
xlabel(ax, 'Y  (mm)');
ylabel(ax, 'Z (mm)');
if isempty(ratioLabel)
    titleStr = sprintf('Mid-plane YZ curvature — %s', thetaLabel);
else
    titleStr = sprintf('Mid-plane YZ curvature — %s | %s', ratioLabel, thetaLabel);
end
title(ax, titleStr, 'FontWeight', 'bold', 'Interpreter', 'tex');
legend(ax, 'Location', 'best');

if ismember('X', solid.Properties.VariableNames)
    xVal = mean(double(solid.X), 'omitnan');
    subtitle(ax, sprintf('X = %.4g mm (%d nodes)', xVal, numel(Y0s)));
else
    subtitle(ax, sprintf('%d nodes', numel(Y0s)));
end

exportgraphics(fig, fullfile(outDir, 'curvature_yz_overlay.png'), 'Resolution', 300);
close(fig);

end

function out = format_theta_label(thetaDeg)
if isnan(thetaDeg)
    out = '\theta = N/A';
else
    out = sprintf('\\theta = %d°', round(thetaDeg));
end
end

function [yArc, zArc] = build_theory_arc(Y0s, Z0s, Ysol, Ysh, Zsol, theory)
yArc = [];
zArc = [];

if ~isfield(theory, 'kappa') || isnan(theory.kappa) || theory.kappa <= 0
    return;
end

R = 1 / theory.kappa;
yDefAll = [Ysol(:); Ysh(:)];
yMin = min(yDefAll, [], 'omitnan');
yMax = max(yDefAll, [], 'omitnan');
Delta = yMax - yMin;
if ~(Delta > 0) || ~(R > Delta / 2)
    return;
end

yMid = 0.5 * (yMin + yMax);
zChord = mean([Z0s(1), Z0s(end)], 'omitnan');
delta = R - sqrt(max(0, R^2 - (Delta / 2)^2));

% Match arc orientation to measured displacement sign.
dzMean = mean(Zsol - Z0s, 'omitnan');

%  Ta negativt her for å flippe kurven 
signDir = -sign(dzMean);
if signDir == 0
    signDir = -1;
end

zCenter = zChord - signDir * (R - delta);
yArc = linspace(yMin, yMax, 400).';
radTerm = max(0, R^2 - (yArc - yMid).^2);
zArc = zCenter + signDir * sqrt(radTerm);
end
