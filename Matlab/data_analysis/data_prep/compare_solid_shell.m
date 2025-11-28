function [metrics, errorTable] = compare_solid_shell(solid, shell)
%COMPARE_SOLID_SHELL Compute node-wise errors and aggregate metrics.

arguments
    solid table
    shell table
end

fields = {'U1','U2','U3','S11','S22','S33','S12','S13','S23','SMises'};

% Ensure base coordinate columns are numeric.
baseCols = {'Label','X','Y','Z'};
for i = 1:numel(baseCols)
    c = baseCols{i};
    if ~isnumeric(solid.(c))
        error('compare_solid_shell:NonNumeric', 'Column %s in solid is not numeric (%s)', c, class(solid.(c)));
    end
    if ~isnumeric(shell.(c))
        error('compare_solid_shell:NonNumeric', 'Column %s in shell is not numeric (%s)', c, class(shell.(c)));
    end
end

metrics = struct();
errorTable = table(double(solid.Label), double(solid.X), double(solid.Y), double(solid.Z), ...
    'VariableNames', {'Label', 'X', 'Y', 'Z'});

for iField = 1:numel(fields)
    name = fields{iField};
    if ~ismember(name, solid.Properties.VariableNames) || ~ismember(name, shell.Properties.VariableNames)
        error('compare_solid_shell:MissingField', 'Field %s not present in both tables.', name);
    end

    sVal = double(solid.(name));
    shVal = double(shell.(name));

    err = shVal - sVal;
    scale = max(abs(sVal));
    relEps = max(1e-12, 1e-6 * scale);
    relErr = err ./ (abs(sVal) + relEps);

    mae = mean(abs(err), 'omitnan');
    rmse = sqrt(mean(err.^2, 'omitnan'));
    maxAbs = max(abs(err), [], 'omitnan');
    meanRel = mean(abs(relErr), 'omitnan');
    maxRel = max(abs(relErr), [], 'omitnan');

    R = NaN;
    solidVar = sVal;
    shellVar = shVal;
    if std(solidVar) > 0 && std(shellVar) > 0
        cc = corrcoef(solidVar, shellVar);
        if size(cc, 1) >= 2
            R = cc(1, 2);
        end
    end

    metrics.(name) = struct( ...
        'MAE', mae, ...
        'RMSE', rmse, ...
        'MaxAbs', maxAbs, ...
        'MeanRel', meanRel, ...
        'MaxRel', maxRel, ...
        'R', R);

    errorTable.(sprintf('d_%s', name)) = err;
    errorTable.(sprintf('rel_%s', name)) = relErr;
end

end
