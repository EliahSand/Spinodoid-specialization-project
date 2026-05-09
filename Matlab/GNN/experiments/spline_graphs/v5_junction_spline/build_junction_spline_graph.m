function jlGraph = build_junction_spline_graph(skelMask, radiusMap, boundaryMask, xVals, yVals, varargin)
%BUILD_JUNCTION_SPLINE_GRAPH  Junction-spline (v5) graph from skeleton mask.
%
%   Nodes  = medial-axis junctions (deg>=3) and endpoints (deg==1).
%   Edges  = one per ligament (multi-edges preserved, self-loops preserved).
%   Node features (4): [x, y, radius, is_boundary]
%   Edge features (4): [phys_length, mean_radius, cos2theta, sin2theta]
%
%   The arc-length and arc-weighted mean radius are computed analytically
%   from a cubic B-spline fit to the ligament pixel polyline, which filters
%   8-connected staircase noise. Orientation uses the chord direction in a
%   pi-periodic encoding consistent with the global ang_deg convention.
%
%   Inputs:
%     skelMask    - H x W logical, skeleton pixels
%     radiusMap   - H x W double, physical radius at each pixel
%     boundaryMask- H x W logical, domain-boundary flag per grid cell
%     xVals       - 1 x W physical x for each column
%     yVals       - 1 x H physical y for each row
%   Optional name-value:
%     ErrorTolerance           (default 2*grid_spacing)
%     MaxControlPointsPerBranch (default 8)
%     MinControlPointsPerBranch (default 2)

p = inputParser;
p.addParameter('ErrorTolerance',            [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
p.addParameter('MaxControlPointsPerBranch', 8,  @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('MinControlPointsPerBranch', 2,  @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.parse(varargin{:});
opts = p.Results;

skelMask     = logical(skelMask);
radiusMap    = double(radiusMap);
if nargin >= 3 && ~isempty(boundaryMask)
    boundaryMask = logical(boundaryMask);
else
    boundaryMask = false(size(skelMask));
end
xVals        = double(xVals(:));
yVals        = double(yVals(:));

[H, W] = size(skelMask);
h = estimate_h(xVals, yVals);
if isempty(opts.ErrorTolerance)
    opts.ErrorTolerance = 2 * h;
end

% -------------------------------------------------------------------------
% 1. Physical coordinates + radii per skeleton pixel
% -------------------------------------------------------------------------
skelLin = find(skelMask);
nSkel   = numel(skelLin);

if nSkel == 0
    jlGraph = empty_graph();
    return;
end

localAt = zeros(H, W);
localAt(skelLin) = 1:nSkel;

[rSkel, cSkel] = ind2sub([H, W], skelLin);
xSkel  = xVals(cSkel);
ySkel  = yVals(rSkel);
rhoSkel = radiusMap(skelLin);
bndSkel = boundaryMask(skelLin);

% -------------------------------------------------------------------------
% 2. 8-connected degree + undirected edge list
% -------------------------------------------------------------------------
drs = [-1 -1 -1  0  0  1  1  1];
dcs = [-1  0  1 -1  1 -1  0  1];

deg     = zeros(nSkel, 1);
nEdgeBuf = 4 * nSkel + 8;
ePairs   = zeros(nEdgeBuf, 2);
nEdge    = 0;
for i = 1:nSkel
    r = rSkel(i);
    c = cSkel(i);
    for d = 1:8
        rr = r + drs(d);
        cc = c + dcs(d);
        if rr < 1 || rr > H || cc < 1 || cc > W, continue; end
        j = localAt(rr, cc);
        if j == 0, continue; end
        deg(i) = deg(i) + 1;
        if j > i
            nEdge = nEdge + 1;
            if nEdge > nEdgeBuf
                ePairs   = [ePairs; zeros(nEdgeBuf, 2)]; %#ok<AGROW>
                nEdgeBuf = nEdgeBuf * 2;
            end
            ePairs(nEdge, :) = [i, j];
        end
    end
end
ePairs = ePairs(1:nEdge, :);

% -------------------------------------------------------------------------
% 3. Per-node neighbour + edge lists (O(nEdge) build)
% -------------------------------------------------------------------------
nbrs  = cell(nSkel, 1);
nedge = cell(nSkel, 1);
for e = 1:nEdge
    u = ePairs(e,1); v = ePairs(e,2);
    nbrs{u}(end+1) = v;  nedge{u}(end+1) = e;
    nbrs{v}(end+1) = u;  nedge{v}(end+1) = e;
end

% -------------------------------------------------------------------------
% 4. Classify pixels
% -------------------------------------------------------------------------
isKey     = deg >= 3 | deg == 1 | deg == 0;
keyList   = find(isKey);
nJL       = numel(keyList);
skelToJL  = zeros(nSkel, 1);
skelToJL(keyList) = 1:nJL;

% -------------------------------------------------------------------------
% 5. Trace ligaments
% -------------------------------------------------------------------------
visitedEdge = false(nEdge, 1);

maxLig = nEdge + nJL + 1;
jl_ei  = zeros(maxLig, 2);
jl_ef  = zeros(maxLig, 4);
nJLEdge = 0;

for ki = 1:nJL
    u = keyList(ki);
    for ni = 1:numel(nbrs{u})
        v = nbrs{u}(ni);
        e = nedge{u}(ni);
        if visitedEdge(e), continue; end
        path = trace_lig(u, v, nbrs, nedge, deg, visitedEdge);
        eIds = path_eids(path, nbrs, nedge);
        visitedEdge(eIds) = true;
        srcJL = skelToJL(path(1));
        dstJL = skelToJL(path(end));
        if srcJL == 0 || dstJL == 0, continue; end
        feat = edge_features(path, xSkel, ySkel, rhoSkel, opts);
        nJLEdge = nJLEdge + 1;
        jl_ei(nJLEdge,:) = [srcJL, dstJL];
        jl_ef(nJLEdge,:) = feat;
    end
end

% Isolated pixels (deg==0) become self-loop edges with a zero-length feature
for ki = 1:nJL
    u = keyList(ki);
    if deg(u) == 0
        nJLEdge = nJLEdge + 1;
        jlIdx = skelToJL(u);
        jl_ei(nJLEdge,:) = [jlIdx, jlIdx];
        jl_ef(nJLEdge,:) = [0, rhoSkel(u), 0, 0];
    end
end

% Pure cycles: unvisited edges belong to all-degree-2 components
unvis = find(~visitedEdge);
while ~isempty(unvis)
    u = ePairs(unvis(1), 1);
    % Add a synthetic JL node for this cycle's anchor
    nJL = nJL + 1;
    skelToJL(u) = nJL;
    keyList(nJL) = u;    % u is already a valid skeleton pixel index
    % Trace the cycle
    startN = nbrs{u}(1);
    path = trace_lig(u, startN, nbrs, nedge, deg, visitedEdge);
    eIds = path_eids(path, nbrs, nedge);
    visitedEdge(eIds) = true;
    feat = edge_features(path, xSkel, ySkel, rhoSkel, opts);
    nJLEdge = nJLEdge + 1;
    jl_ei(nJLEdge,:) = [nJL, nJL];
    jl_ef(nJLEdge,:) = feat;
    unvis = find(~visitedEdge);
end

jl_ei = jl_ei(1:nJLEdge, :);
jl_ef = jl_ef(1:nJLEdge, :);

% -------------------------------------------------------------------------
% 6. Build output
% -------------------------------------------------------------------------
jlCoords  = [xSkel(keyList), ySkel(keyList)];
jlRadius  = rhoSkel(keyList);
jlIsBnd   = double(bndSkel(keyList));

node_features = [jlCoords, jlRadius(:), jlIsBnd(:)];  % nJL x 4

jlGraph = struct();
jlGraph.num_nodes        = nJL;
jlGraph.num_edges        = nJLEdge;
jlGraph.node_coords      = jlCoords;
jlGraph.node_features    = node_features;
jlGraph.feature_names    = {'x', 'y', 'radius', 'is_boundary'};
jlGraph.edges_local      = jl_ei;
jlGraph.edge_index       = int32(jl_ei.');
jlGraph.edge_features    = jl_ef;
jlGraph.edge_feature_names = {'phys_length', 'mean_radius', 'cos2theta', 'sin2theta'};
jlGraph.schema_version   = 5;
jlGraph.kind             = 'junction_spline_graph';
end

% =========================================================================
% Helpers
% =========================================================================

function path = trace_lig(startNode, nextNode, nbrs, nedge, deg, visitedEdge)
path = [startNode, nextNode];
prev = startNode;
curr = nextNode;
while true
    if curr == startNode && numel(path) > 2, break; end  % closed cycle
    if deg(curr) ~= 2, break; end                         % hit key node
    ns = nbrs{curr};
    es = nedge{curr};
    found = false;
    for k = 1:numel(ns)
        nxt = ns(k);
        e   = es(k);
        if nxt == prev, continue; end
        if visitedEdge(e) && nxt ~= startNode, continue; end
        path(end+1) = nxt; %#ok<AGROW>
        prev = curr;
        curr = nxt;
        found = true;
        break;
    end
    if ~found, break; end
end
end

function eIds = path_eids(path, nbrs, nedge)
eIds = zeros(numel(path)-1, 1);
for k = 1:numel(path)-1
    u = path(k); v = path(k+1);
    ns = nbrs{u}; es = nedge{u};
    for m = 1:numel(ns)
        if ns(m) == v
            eIds(k) = es(m);
            break;
        end
    end
end
eIds = eIds(eIds > 0);
end

function feat = edge_features(path, xSkel, ySkel, rhoSkel, opts)
pts3 = [xSkel(path(:)), ySkel(path(:)), rhoSkel(path(:))];

% Chord orientation (pi-periodic)
dx = pts3(end,1) - pts3(1,1);
dy = pts3(end,2) - pts3(1,2);
theta = atan2(dy, dx);
cos2t = cos(2 * theta);
sin2t = sin(2 * theta);

if size(pts3, 1) < 2
    feat = [0, rhoSkel(path(1)), cos2t, sin2t];
    return;
end

% Fit cubic B-spline for smoothed length + mean radius
fit = fit_branch_spline(pts3, opts.ErrorTolerance, ...
    opts.MinControlPointsPerBranch, opts.MaxControlPointsPerBranch);
[L, meanR] = spline_aggregates(fit.control_points, fit.knots, fit.degree);

feat = [L, meanR, cos2t, sin2t];
end

function [L, meanR] = spline_aggregates(C, knots, degree)
% 8-point Gauss-Legendre quadrature on [0,1]
tGL = [0.019855071751232, 0.101666776969940, 0.237233795041836, 0.408282678752175, ...
       0.591717321247825, 0.762766204958164, 0.898333223030060, 0.980144928248768];
wGL = [0.050614268145185, 0.111190517226687, 0.156853322938944, 0.181341891689181, ...
       0.181341891689181, 0.156853322938944, 0.111190517226687, 0.050614268145185];

nCtrl = size(C, 1);

% Build basis + derivative basis at GL nodes
B  = zeros(8, nCtrl);
dB = zeros(8, nCtrl);
for i = 1:nCtrl
    [Bi, dBi] = bspline_basis_and_deriv(i, degree, tGL(:), knots, nCtrl);
    B(:,i)  = Bi;
    dB(:,i) = dBi;
end
B(tGL == 1, :) = 0; B(tGL == 1, nCtrl) = 1;

Ceval  = B  * C;   % 8 x 3  (x,y,r at GL nodes)
dCeval = dB * C;   % 8 x 3  (dx/dt, dy/dt, dr/dt)

speeds = sqrt(dCeval(:,1).^2 + dCeval(:,2).^2);  % 8 x 1

L     = dot(wGL(:), speeds);
if L < eps
    meanR = mean(C(:,3));
    return;
end
meanR = dot(wGL(:), speeds .* max(0, Ceval(:,3))) / L;
end

function [N, dN] = bspline_basis_and_deriv(i, degree, t, knots, nCtrl)
N  = bspline_basis(i, degree, t, knots, nCtrl);
if degree == 0
    dN = zeros(size(t));
    return;
end
denom1 = knots(i + degree)     - knots(i);
denom2 = knots(i + degree + 1) - knots(i + 1);
dN = zeros(size(t));
if denom1 > eps
    dN = dN + degree / denom1 * bspline_basis(i,   degree-1, t, knots, nCtrl);
end
if denom2 > eps
    dN = dN - degree / denom2 * bspline_basis(i+1, degree-1, t, knots, nCtrl);
end
end

% ---- Spline fitting (self-contained copy so v4 code is untouched) --------

function fit = fit_branch_spline(pts3, errorTol, minCtrl, maxCtrl)
nData = size(pts3, 1);
if nData <= 1
    fit = struct('control_points', pts3, 'knots', [0 1], 'degree', 0, 'fit_error', 0);
    return;
end
t = chord_params(pts3(:,1:2));
maxCtrl = min(maxCtrl, nData);
minCtrl = min(max(2, minCtrl), maxCtrl);
best = [];
for nCtrl = minCtrl:maxCtrl
    deg  = min(3, nCtrl - 1);
    kv   = open_uniform_knots(nCtrl, deg);
    B    = bspline_basis_matrix(t, nCtrl, deg, kv);
    Cfit = fit_control_points(B, pts3);
    approx = B * Cfit;
    err  = max(sqrt(sum((pts3 - approx).^2, 2)));
    best = struct('control_points', Cfit, 'knots', kv, 'degree', deg, 'fit_error', err);
    if err <= errorTol + eps, break; end
end
fit = best;
end

function t = chord_params(xy)
if size(xy, 1) == 1, t = 0; return; end
ds = sqrt(sum(diff(xy,1,1).^2, 2));
s  = [0; cumsum(ds)];
if s(end) <= eps
    t = linspace(0,1,size(xy,1)).';
else
    t = s / s(end);
end
t(1) = 0; t(end) = 1;
end

function kv = open_uniform_knots(nCtrl, degree)
ni = nCtrl - degree - 1;
if ni > 0; interior = (1:ni)/(ni+1); else; interior = zeros(1,0); end
kv = [zeros(1,degree+1), interior, ones(1,degree+1)];
end

function B = bspline_basis_matrix(t, nCtrl, degree, kv)
t = t(:);
B = zeros(numel(t), nCtrl);
for i = 1:nCtrl
    B(:,i) = bspline_basis(i, degree, t, kv, nCtrl);
end
B(t == 1, :) = 0; B(t == 1, nCtrl) = 1;
end

function N = bspline_basis(i, degree, t, kv, nCtrl)
if degree == 0
    N = double(kv(i) <= t & t < kv(i+1));
    if i == nCtrl, N(t == 1) = 1; end
    return;
end
da = kv(i+degree)   - kv(i);
db = kv(i+degree+1) - kv(i+1);
A = zeros(size(t)); B = zeros(size(t));
if da > eps, A = ((t - kv(i))   / da) .* bspline_basis(i,   degree-1, t, kv, nCtrl); end
if db > eps, B = ((kv(i+degree+1) - t) / db) .* bspline_basis(i+1, degree-1, t, kv, nCtrl); end
N = A + B;
end

function C = fit_control_points(B, pts3)
nCtrl = size(B,2);
if nCtrl <= 2
    C = B \ pts3;
    return;
end
C = zeros(nCtrl, size(pts3,2));
C(1,:)   = pts3(1,:);
C(end,:) = pts3(end,:);
rhs = pts3 - B(:,[1,nCtrl]) * C([1,nCtrl],:);
C(2:end-1,:) = B(:,2:end-1) \ rhs;
end

function h = estimate_h(xVals, yVals)
dx = diff(sort(xVals(:))); dy = diff(sort(yVals(:)));
vals = [dx(dx>0); dy(dy>0)];
if isempty(vals), h = 1; else, h = median(vals); end
end

function g = empty_graph()
g = struct('num_nodes', 0, 'num_edges', 0, ...
    'node_coords', zeros(0,2), 'node_features', zeros(0,4), ...
    'feature_names', {{'x','y','radius','is_boundary'}}, ...
    'edges_local', zeros(0,2), 'edge_index', int32(zeros(2,0)), ...
    'edge_features', zeros(0,4), ...
    'edge_feature_names', {{'phys_length','mean_radius','cos2theta','sin2theta'}}, ...
    'schema_version', 5, 'kind', 'junction_spline_graph');
end
