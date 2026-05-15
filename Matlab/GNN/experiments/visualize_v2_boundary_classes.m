% visualize_v2_boundary_classes
%
% Quick sanity-check for the V2 boundary one-hot. Loads one sample.mat from
% data/dataset_hybrid_d100_v2/samples and scatter-plots nodes coloured by
% their class index {interior, clamped_left, loaded_right, free_top_bot}.
%
% Edit SAMPLE_DIR / SAMPLE_NAME (or leave SAMPLE_NAME = '' to auto-pick the
% first available sample) and run.

scriptDir = fileparts(mfilename('fullpath'));
gnnRoot   = fileparts(scriptDir);
SAMPLES_DIR = fullfile(gnnRoot, 'data', 'dataset_hybrid_d100_v2', 'samples');
SAMPLE_NAME = '';   % e.g. 'tr50_ang030_lamellar_N128_1x1'; '' picks first

if isempty(SAMPLE_NAME)
    items = dir(SAMPLES_DIR);
    items = items([items.isdir] & ~startsWith({items.name}, '.'));
    if isempty(items)
        error('No sample folders found under %s', SAMPLES_DIR);
    end
    % Walk into the first folder that contains a sample.mat (handles nested layouts).
    SAMPLE_NAME = items(1).name;
    while ~isfile(fullfile(SAMPLES_DIR, SAMPLE_NAME, 'sample.mat'))
        sub = dir(fullfile(SAMPLES_DIR, SAMPLE_NAME));
        sub = sub([sub.isdir] & ~startsWith({sub.name}, '.'));
        if isempty(sub), error('No sample.mat under %s', fullfile(SAMPLES_DIR, SAMPLE_NAME)); end
        SAMPLE_NAME = fullfile(SAMPLE_NAME, sub(1).name);
    end
end

matPath = fullfile(SAMPLES_DIR, SAMPLE_NAME, 'sample.mat');
fprintf('Loading: %s\n', matPath);
d = load(matPath);
gd = d.gnn_data;

assert(strcmpi(string(gd.schema_version_label), 'V2'), ...
    'Expected V2 sample, got schema_version_label=%s', string(gd.schema_version_label));
assert(size(gd.x, 2) == 7, 'Expected 7-dim node features; got %d.', size(gd.x, 2));

xy = gd.x(:, 1:2);
oneHot = gd.x(:, 4:7);

[~, classIdx] = max(oneHot, [], 2);   % 1=interior 2=clamped 3=loaded 4=free

classNames = {'interior', 'clamped\_left', 'loaded\_right', 'free\_top\_bot'};
classCols  = [0.75 0.75 0.75;
              0.85 0.10 0.10;
              0.10 0.45 0.85;
              0.20 0.65 0.20];

figure('Color', 'w', 'Position', [100 100 720 640]); hold on;
for c = 1:4
    sel = classIdx == c;
    scatter(xy(sel, 1), xy(sel, 2), 22, classCols(c, :), 'filled', ...
        'DisplayName', sprintf('%s (n=%d)', classNames{c}, nnz(sel)));
end
axis equal tight; box on;
xlabel('x'); ylabel('y');
title(sprintf('V2 boundary one-hot — %s', strrep(SAMPLE_NAME, '_', '\_')));
legend('Location', 'eastoutside');

fprintf('Class counts: interior=%d  clamped_left=%d  loaded_right=%d  free_top_bot=%d  (total=%d)\n', ...
    nnz(classIdx==1), nnz(classIdx==2), nnz(classIdx==3), nnz(classIdx==4), numel(classIdx));
