function pathOut = unique_path(pathIn)
pathOut = pathIn;
if isfile(pathOut)
    [folder, base, ext] = fileparts(pathOut);
    v = 1;
    while isfile(fullfile(folder, sprintf('%s_v%d%s', base, v, ext)))
        v = v + 1;
    end
    pathOut = fullfile(folder, sprintf('%s_v%d%s', base, v, ext));
end
end
