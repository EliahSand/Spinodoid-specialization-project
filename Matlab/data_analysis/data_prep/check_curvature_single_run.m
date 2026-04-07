function [report, reportTable] = check_curvature_single_run(runFolder, solid, shell, opts)
%CHECK_CURVATURE_SINGLE_RUN Curvature check for a single solid/shell run.
%   Implements the laminate-theory curvature relation used in
%   LaminateTheory/curveCalculation.ipynb and compares it with curvature
%   measured from the deformed YZ profile of solid and shell outputs.
%
%   [report, reportTable] = check_curvature_single_run(runFolder, solid, shell)
%
%   Outputs:
%     report      - struct with theory / solid / shell curvature details.
%     reportTable - tidy table for CSV export.

arguments
    runFolder (1, :) char
    solid table
    shell table
    opts.Poisson (1, 1) double = 0.4
    opts.DispFraction (1, 1) double = 0.1
end

runLogPath = fullfile(runFolder, 'run_log.txt');
logData = parse_run_log(runLogPath);

theory = struct( ...
    't', NaN, ...
    'h', NaN, ...
    'w', NaN, ...
    'd', NaN, ...
    'phi', NaN, ...
    'nu', opts.Poisson, ...
    'dispFraction', opts.DispFraction, ...
    'kappa', NaN, ...
    'radius', NaN, ...
    'hasInputs', false);

if isfield(logData, 't_base_m') && isfield(logData, 't_spin_m') ...
        && isfield(logData, 'line_width_m') && isfield(logData, 'line_gap_factor')
    theory.t = logData.t_base_m;
    theory.h = logData.t_spin_m;
    theory.w = logData.line_width_m;
    theory.d = theory.w * logData.line_gap_factor;
    if theory.w > 0 && theory.d > 0
        theory.phi = theory.w / (theory.w + theory.d);
        num = 6 * theory.t * theory.h * theory.phi * (theory.t + theory.h);
        den = (theory.t^4) + (4 * theory.t^3 * theory.h * theory.phi) ...
            + (6 * theory.t^2 * theory.h^2 * theory.phi) ...
            + (4 * theory.t * theory.h^3 * theory.phi) ...
            + (theory.h^4 * theory.phi^2);
        if den > 0
            theory.kappa = (num / den) * theory.nu * theory.dispFraction;
            if theory.kappa > 0
                theory.radius = 1 / theory.kappa;
            end
            theory.hasInputs = true;
        end
    end
end

solidMeas = measure_profile_curvature(solid);
shellMeas = measure_profile_curvature(shell);

epsDen = max(1e-12, abs(theory.kappa));
solidRelTheory = NaN;
shellRelTheory = NaN;
if ~isnan(theory.kappa)
    solidRelTheory = abs(solidMeas.kappaCircle - theory.kappa) / epsDen;
    shellRelTheory = abs(shellMeas.kappaCircle - theory.kappa) / epsDen;
end

report = struct();
report.runFolder = runFolder;
report.runLogPath = runLogPath;
report.theory = theory;
report.solid = solidMeas;
report.shell = shellMeas;
report.solid.RelErrTheory = solidRelTheory;
report.shell.RelErrTheory = shellRelTheory;

models = ["theory"; "solid"; "shell"];
kappa = [theory.kappa; solidMeas.kappaCircle; shellMeas.kappaCircle];
radius = [theory.radius; solidMeas.radiusCircle; shellMeas.radiusCircle];
delta = [NaN; solidMeas.delta; shellMeas.delta];
chord = [NaN; solidMeas.Delta; shellMeas.Delta];
kappaSag = [NaN; solidMeas.kappaSagitta; shellMeas.kappaSagitta];
fitRmse = [NaN; solidMeas.circleFitRMSE; shellMeas.circleFitRMSE];
relTheory = [NaN; solidRelTheory; shellRelTheory];

reportTable = table(models, kappa, radius, delta, chord, kappaSag, fitRmse, relTheory, ...
    'VariableNames', {'Model', 'Kappa_1_per_m', 'Radius_m', 'Sagitta_m', ...
    'Chord_m', 'KappaSagitta_1_per_m', 'CircleFitRMSE_m', 'RelErrToTheory'});

end

function out = measure_profile_curvature(tbl)
Y = double(tbl.Y) + double(tbl.U2);
Z = double(tbl.Z) + double(tbl.U3);

valid = isfinite(Y) & isfinite(Z);
Y = Y(valid);
Z = Z(valid);

out = struct( ...
    'nPoints', numel(Y), ...
    'Delta', NaN, ...
    'delta', NaN, ...
    'kappaSagitta', NaN, ...
    'radiusSagitta', NaN, ...
    'kappaCircle', NaN, ...
    'radiusCircle', NaN, ...
    'circleCenterY', NaN, ...
    'circleCenterZ', NaN, ...
    'circleFitRMSE', NaN);

if numel(Y) < 3
    return;
end

[Y, idx] = sort(Y);
Z = Z(idx);

out.Delta = max(Y) - min(Y);
if out.Delta <= 0
    return;
end

% Sagitta from max perpendicular distance to the chord.
p1 = [Y(1), Z(1)];
p2 = [Y(end), Z(end)];
v = p2 - p1;
vn = hypot(v(1), v(2));
if vn > 0
    signedDist = ((Y - p1(1)) * v(2) - (Z - p1(2)) * v(1)) / vn;
    out.delta = max(abs(signedDist), [], 'omitnan');
end

if ~isnan(out.delta) && out.delta > 0
    out.radiusSagitta = (out.Delta^2) / (8 * out.delta) + out.delta / 2;
    if out.radiusSagitta > 0
        out.kappaSagitta = 1 / out.radiusSagitta;
    end
end

% Algebraic circle fit (Kasa).
A = [2 * Y, 2 * Z, ones(size(Y))];
b = Y.^2 + Z.^2;
x = A \ b;
yc = x(1);
zc = x(2);
c0 = x(3);
r2 = c0 + yc^2 + zc^2;
if isfinite(r2) && r2 > 0
    R = sqrt(r2);
    rr = hypot(Y - yc, Z - zc) - R;
    out.radiusCircle = R;
    out.kappaCircle = 1 / R;
    out.circleCenterY = yc;
    out.circleCenterZ = zc;
    out.circleFitRMSE = sqrt(mean(rr.^2, 'omitnan'));
end

end

function s = parse_run_log(runLogPath)
s = struct();
if ~isfile(runLogPath)
    return;
end

fid = fopen(runLogPath, 'r');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    parts = split(string(line), ":");
    if numel(parts) < 2
        continue;
    end
    key = strtrim(lower(parts(1)));
    valStr = strtrim(join(parts(2:end), ":"));
    val = str2double(valStr);
    if isnan(val)
        continue;
    end
    switch key
        case "t_base_m"
            s.t_base_m = val;
        case "t_spin_m"
            s.t_spin_m = val;
        case "line_width_m"
            s.line_width_m = val;
        case "line_gap_factor"
            s.line_gap_factor = val;
    end
end

end
