function metricsTable = metrics_to_table(metrics)
%METRICS_TO_TABLE Convert metrics struct into a tidy table for export.

fields = fieldnames(metrics);
n = numel(fields);
metricsTable = table('Size', [n, 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'Field','MAE','RMSE','MaxAbs','MeanRel','MaxRel','R'});

for i = 1:n
    fname = fields{i};
    m = metrics.(fname);
    metricsTable.Field(i) = string(fname);
    metricsTable.MAE(i) = double(m.MAE);
    metricsTable.RMSE(i) = double(m.RMSE);
    metricsTable.MaxAbs(i) = double(m.MaxAbs);
    metricsTable.MeanRel(i) = double(m.MeanRel);
    metricsTable.MaxRel(i) = double(m.MaxRel);
    metricsTable.R(i) = double(m.R);
end
end
