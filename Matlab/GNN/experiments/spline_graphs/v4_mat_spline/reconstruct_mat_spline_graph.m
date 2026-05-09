function [reconMask, sampled, metrics] = reconstruct_mat_spline_graph(matGraph, xVals, yVals, varargin)
%RECONSTRUCT_MAT_SPLINE_GRAPH Rebuild a mask from MAT/spline branch circles.
%   [reconMask, sampled, metrics] = RECONSTRUCT_MAT_SPLINE_GRAPH(matGraph, xVals, yVals)
%   densely samples each stored branch spline as (x,y,radius) and unions the
%   corresponding medial disks on the native grid.

    p = inputParser;
    p.addParameter('SamplesPerGridStep', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('MinSamplesPerBranch', 32, @(x) isnumeric(x) && isscalar(x) && x >= 2);
    p.addParameter('OccupancyMask', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
    p.parse(varargin{:});
    opts = p.Results;

    xVals = double(xVals(:));
    yVals = double(yVals(:));
    if isempty(xVals) || isempty(yVals)
        error('reconstruct_mat_spline_graph:BadGrid', 'xVals and yVals must be nonempty.');
    end

    nodeFeatures = get_node_features(matGraph);
    if size(nodeFeatures, 2) < 3
        error('reconstruct_mat_spline_graph:BadFeatures', ...
            'matGraph node features must contain at least [x,y,radius].');
    end
    if ~isfield(matGraph, 'branches')
        error('reconstruct_mat_spline_graph:MissingBranches', ...
            'matGraph is missing branch metadata.');
    end

    h = estimate_grid_spacing(xVals, yVals);
    branches = matGraph.branches;
    points = zeros(0, 3);
    branchIndex = zeros(0, 1);

    for b = 1:numel(branches)
        if ~isfield(branches(b), 'control_indices')
            continue;
        end
        ctrlIds = double(branches(b).control_indices(:));
        ctrlIds = ctrlIds(ctrlIds >= 1 & ctrlIds <= size(nodeFeatures, 1));
        if isempty(ctrlIds)
            continue;
        end

        ctrl = double(nodeFeatures(ctrlIds, 1:3));
        ctrl(:, 3) = max(0, ctrl(:, 3));
        isCycle = isfield(branches(b), 'is_cycle') && logical(branches(b).is_cycle);
        len = branch_length_guess(ctrl(:, 1:2), branches(b), isCycle);
        nSamples = max([opts.MinSamplesPerBranch, ceil(len / max(h, eps) * opts.SamplesPerGridStep), 8 * size(ctrl, 1)]);

        curve = evaluate_branch_curve_3d(ctrl, branches(b), nSamples);
        if isempty(curve)
            curve = sample_control_polyline(ctrl, nSamples, isCycle);
        end
        if isempty(curve)
            continue;
        end

        points = [points; curve]; %#ok<AGROW>
        branchIndex = [branchIndex; repmat(b, size(curve, 1), 1)]; %#ok<AGROW>
    end

    reconMask = rasterize_medial_disks(points, xVals, yVals);

    sampled = struct();
    sampled.points = points;
    sampled.branch_index = branchIndex;
    sampled.grid_spacing = h;
    sampled.sample_count = size(points, 1);
    sampled.samples_per_grid_step = opts.SamplesPerGridStep;

    metrics = struct();
    if ~isempty(opts.OccupancyMask)
        occMask = logical(opts.OccupancyMask);
        if ~isequal(size(occMask), size(reconMask))
            error('reconstruct_mat_spline_graph:BadOccupancySize', ...
                'OccupancyMask must be size numel(yVals) x numel(xVals).');
        end
        metrics = reconstruction_metrics(occMask, reconMask, h);
        sampled.metrics = metrics;
    end
end

function X = get_node_features(matGraph)
    if isfield(matGraph, 'node_features')
        X = double(matGraph.node_features);
    elseif isfield(matGraph, 'x')
        X = double(matGraph.x);
    else
        error('reconstruct_mat_spline_graph:MissingFeatures', ...
            'matGraph must contain node_features or x.');
    end

    n = get_field_default(matGraph, 'num_nodes', NaN);
    if ~isnan(n) && size(X, 2) == n && size(X, 1) ~= n
        X = X.';
    end
end

function h = estimate_grid_spacing(xVals, yVals)
    dx = diff(sort(xVals(:)));
    dy = diff(sort(yVals(:)));
    vals = [dx(dx > 0); dy(dy > 0)];
    if isempty(vals)
        h = 1;
    else
        h = median(vals);
    end
end

function len = branch_length_guess(xy, branch, isCycle)
    if isfield(branch, 'path_length') && isfinite(branch.path_length) && branch.path_length > 0
        len = double(branch.path_length);
        return;
    end
    if size(xy, 1) < 2
        len = 0;
        return;
    end
    if isCycle
        xy = [xy; xy(1, :)]; %#ok<AGROW>
    end
    len = sum(sqrt(sum(diff(xy, 1, 1) .^ 2, 2)));
end

function curve = evaluate_branch_curve_3d(ctrl, branch, nSamples)
    curve = [];
    if size(ctrl, 1) < 2 || ~isfield(branch, 'degree') || ~isfield(branch, 'knots')
        return;
    end

    nCtrl = size(ctrl, 1);
    degree = double(branch.degree);
    knots = double(branch.knots(:).');
    if degree < 1 || nCtrl <= degree || numel(knots) ~= nCtrl + degree + 1
        return;
    end

    t = linspace(0, 1, nSamples).';
    B = bspline_basis_matrix(t, nCtrl, degree, knots);
    curve = B * ctrl;
    curve(:, 3) = max(0, curve(:, 3));
end

function B = bspline_basis_matrix(t, nCtrl, degree, knots)
    t = t(:);
    B = zeros(numel(t), nCtrl);
    for i = 1:nCtrl
        B(:, i) = bspline_basis(i, degree, t, knots, nCtrl);
    end
    B(t == 1, :) = 0;
    B(t == 1, nCtrl) = 1;
end

function N = bspline_basis(i, degree, t, knots, nCtrl)
    if degree == 0
        N = double(knots(i) <= t & t < knots(i + 1));
        if i == nCtrl
            N(t == 1) = 1;
        end
        return;
    end

    denomA = knots(i + degree) - knots(i);
    denomB = knots(i + degree + 1) - knots(i + 1);
    A = zeros(size(t));
    B = zeros(size(t));
    if denomA > eps
        A = ((t - knots(i)) / denomA) .* bspline_basis(i, degree - 1, t, knots, nCtrl);
    end
    if denomB > eps
        B = ((knots(i + degree + 1) - t) / denomB) .* bspline_basis(i + 1, degree - 1, t, knots, nCtrl);
    end
    N = A + B;
end

function curve = sample_control_polyline(ctrl, nSamples, isCycle)
    if size(ctrl, 1) == 1
        curve = ctrl;
        return;
    end
    if isCycle
        ctrl = [ctrl; ctrl(1, :)]; %#ok<AGROW>
    end
    ds = sqrt(sum(diff(ctrl(:, 1:2), 1, 1) .^ 2, 2));
    s = [0; cumsum(ds)];
    if s(end) <= eps
        curve = repmat(ctrl(1, :), nSamples, 1);
        return;
    end
    tq = linspace(0, s(end), nSamples).';
    curve = zeros(nSamples, 3);
    for d = 1:3
        curve(:, d) = interp1(s, ctrl(:, d), tq, 'linear');
    end
    curve(:, 3) = max(0, curve(:, 3));
end

function recon = rasterize_medial_disks(points, xVals, yVals)
    recon = false(numel(yVals), numel(xVals));
    if isempty(points)
        return;
    end

    [xGrid, yGrid] = meshgrid(xVals, yVals);
    for i = 1:size(points, 1)
        x = points(i, 1);
        y = points(i, 2);
        r = points(i, 3);
        if ~isfinite(x) || ~isfinite(y) || ~isfinite(r)
            continue;
        end
        if r <= 0
            [~, c] = min(abs(xVals - x));
            [~, row] = min(abs(yVals - y));
            recon(row, c) = true;
            continue;
        end
        cLo = find(xVals >= x - r, 1, 'first');
        cHi = find(xVals <= x + r, 1, 'last');
        rLo = find(yVals >= y - r, 1, 'first');
        rHi = find(yVals <= y + r, 1, 'last');
        if isempty(cLo) || isempty(cHi) || isempty(rLo) || isempty(rHi)
            continue;
        end
        local = (xGrid(rLo:rHi, cLo:cHi) - x) .^ 2 + ...
            (yGrid(rLo:rHi, cLo:cHi) - y) .^ 2 <= r ^ 2 + eps;
        recon(rLo:rHi, cLo:cHi) = recon(rLo:rHi, cLo:cHi) | local;
    end
end

function metrics = reconstruction_metrics(occMask, reconMask, h)
    missing = occMask & ~reconMask;
    extra = reconMask & ~occMask;
    unionMask = occMask | reconMask;
    metrics = struct();
    metrics.missing_pixel_count = nnz(missing);
    metrics.extra_pixel_count = nnz(extra);
    metrics.occupancy_pixel_count = nnz(occMask);
    metrics.reconstruction_pixel_count = nnz(reconMask);
    metrics.iou = nnz(occMask & reconMask) / max(nnz(unionMask), 1);
    metrics.missing_fraction = nnz(missing) / max(nnz(occMask), 1);
    metrics.extra_fraction = nnz(extra) / max(nnz(occMask), 1);

    if ~any(missing(:))
        metrics.one_sided_missing_error = 0;
    elseif ~any(reconMask(:))
        metrics.one_sided_missing_error = inf;
    else
        metrics.one_sided_missing_error = max(bwdist(reconMask, 'euclidean') .* missing, [], 'all') * h;
    end
end

function value = get_field_default(s, fieldName, defaultValue)
    if isstruct(s) && isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end
