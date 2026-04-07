%EXAMPLE_BUILD_GRAPHS Minimal example for full + structural graph generation.

thisDir = fileparts(mfilename('fullpath'));
moduleRoot = fileparts(thisDir);
matlabRoot = fileparts(moduleRoot);

addpath(moduleRoot);
addpath(genpath(fullfile(moduleRoot, 'src')));

inpPath = fullfile(matlabRoot, 'results', 'sheets', 'lamellar', ...
    'sheetCone_tr50_ang000_lamellar_N128_1x1', 'FEA_shell', 'sheet_shell_shell_job.inp');
csvPath = fullfile(matlabRoot, 'results', 'sheets', 'lamellar', ...
    'sheetCone_tr50_ang000_lamellar_N128_1x1', 'FEA_shell', 'midplane_results_shell.csv');

outDir = fullfile(moduleRoot, 'out');
outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'StructuralDetailLevel', 0.25, ...
    'OutDir', outDir, ...
    'Prefix', 'example_spinodal');

disp(outputs.exports);
