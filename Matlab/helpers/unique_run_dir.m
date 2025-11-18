function runDir = unique_run_dir(root, base)
runDir = fullfile(root, base);
suffix = 2;
while exist(runDir, 'dir')
    runDir = fullfile(root, sprintf('%s_run%02d', base, suffix));
    suffix = suffix + 1;
end
if ~exist(runDir, 'dir')
    mkdir(runDir);
end
end
