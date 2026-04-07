function plot_curvature_offset_removed(solid, shell, trLabel, thetaDeg, outDir, opts)
%PLOT_CURVATURE_OFFSET_REMOVED Compare curves using baseline-referenced Z.
%   Plots DeltaZ(Y) = Z_deformed(Y) - Z_baseline(Y) for solid/shell/theory.

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

Y0 = double(solid.Y);
Z0 = double(solid.Z);

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

% Baseline interpolation support.
[Y0s, i0] = sort(Y0);
Z0s = Z0(i0);
[Y0u, iu] = unique(Y0s, 'stable');
Z0u = Z0s(iu);
if numel(Y0u) < 2
    warning('plot_curvature_offset_removed:InsufficientBaseline', ...
        'Not enough distinct baseline Y points to build offset-removed plot.');
    return;
end

[Ysols, is] = sort(Ysol);
Zsols = Zsol(is);
[Yshs, ih] = sort(Ysh);
Zshs = Zsh(ih);

Z0AtSol = interp1(Y0u, Z0u, Ysols, 'linear', 'extrap');
Z0AtSh = interp1(Y0u, Z0u, Yshs, 'linear', 'extrap');
dZsol = Zsols - Z0AtSol;
dZsh = Zshs - Z0AtSh;
[dZsol, ~] = remove_constant_offset(dZsol);
[dZsh, ~] = remove_constant_offset(dZsh);

ratioLabel = format_ratio_label(trLabel);
thetaLabel = format_theta_label(thetaDeg);

fig = figure('Visible', 'off');
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
plot(ax, Ysols, dZsol, '-r', 'LineWidth', 1.5, 'DisplayName', 'Solid');
plot(ax, Yshs, dZsh, '-b', 'LineWidth', 1.5, 'DisplayName', 'Shell');

if isfield(opts.CurvatureCheck, 'theory')
    [yTheory, zTheory] = build_theory_arc(Y0s, Z0s, Ysols, Yshs, Zsols, Zshs, opts.CurvatureCheck.theory);
    if ~isempty(yTheory)
        z0AtTheory = interp1(Y0u, Z0u, yTheory, 'linear', 'extrap');
        dZtheory = zTheory - z0AtTheory;
        [dZtheory, ~] = remove_constant_offset(dZtheory);
        plot(ax, yTheory, dZtheory, '--g', 'LineWidth', 1.5, 'DisplayName', 'Theory');
    end
end

grid(ax, 'on');
xlabel(ax, 'Y');
ylabel(ax, '\DeltaZ (baseline referenced, vertical offset removed)');
title(ax, sprintf('Curvature comparison (baseline + vertical shift removed) — %s | %s', ratioLabel, thetaLabel), ...
    'FontWeight', 'bold', 'Interpreter', 'tex');
legend(ax, 'Location', 'best');

if ismember('X', solid.Properties.VariableNames)
    xVal = mean(double(solid.X), 'omitnan');
    subtitle(ax, sprintf('X = %.4g (%d nodes)', xVal, numel(Y0)));
else
    subtitle(ax, sprintf('%d nodes', numel(Y0)));
end

exportgraphics(fig, fullfile(outDir, 'curvature_yz_offset_removed.png'), 'Resolution', 300);
close(fig);

end

function out = format_theta_label(thetaDeg)
if isnan(thetaDeg)
    out = '\theta = N/A';
else
    out = sprintf('\\theta = %d°', round(thetaDeg));
end
end

function [yArc, zArc] = build_theory_arc(Y0s, Z0s, Ysol, Ysh, Zsol, Zsh, theory)
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

signDir = infer_arc_side(Ysol, Zsol);
if signDir == 0
    signDir = infer_arc_side(Ysh, Zsh);
end
if signDir == 0
    dzMean = mean(Zsol - Z0s, 'omitnan');
    signDir = sign(dzMean);
end
if signDir == 0
    signDir = -1;
end
% Match user convention used in the standard curvature plot.
signDir = -signDir;

zCenter = zChord - signDir * (R - delta);
yArc = linspace(yMin, yMax, 400).';
radTerm = max(0, R^2 - (yArc - yMid).^2);
zArc = zCenter + signDir * sqrt(radTerm);
end

function s = infer_arc_side(Y, Z)
s = 0;
valid = isfinite(Y) & isfinite(Z);
Y = Y(valid);
Z = Z(valid);
if numel(Y) < 3
    return;
end
[Y, idx] = sort(Y);
Z = Z(idx);
p1 = [Y(1), Z(1)];
p2 = [Y(end), Z(end)];
v = p2 - p1;
vn = hypot(v(1), v(2));
if vn <= 0
    return;
end
signedDist = ((Y - p1(1)) * v(2) - (Z - p1(2)) * v(1)) / vn;
s = sign(median(signedDist, 'omitnan'));
end

function [zOut, c] = remove_constant_offset(zIn)
zOut = zIn;
valid = isfinite(zIn);
if ~any(valid)
    c = 0;
    return;
end
c = mean(zIn(valid), 'omitnan');
zOut = zIn - c;
end
