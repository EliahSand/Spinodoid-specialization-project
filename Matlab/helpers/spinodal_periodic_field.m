function [phi, meta] = spinodal_periodic_field(N, L, lambda_vox, bandwidth, nModes, coneDeg)
% Periodic Gaussian random field built from discrete Fourier modes.
% Opposite faces match in value & slope by construction (torus).
% Anisotropy: retain |k| near target shell AND within cones about axes.

% --- Grid in real space (0..L) ------------------------------------------
x = (0:N-1)/N * L;  [X,Y,Z] = ndgrid(x,x,x);

% --- Target |k| in index-units ------------------------------------------
% For a discrete grid, spatial period (in voxels) ~ N/|k_index|.
k_target_idx = N / lambda_vox;                % e.g., lambda=16 -> k≈8 (index units)
k_lo = (1 - bandwidth) * k_target_idx;
k_hi = (1 + bandwidth) * k_target_idx;

% --- Enumerate integer k-vectors on the 3D Fourier lattice --------------
% Note: continuous k = (2π/L) * k_index. Periodicity follows from integer indices.
kk = -floor(N/2):floor((N-1)/2);  % symmetric index set
[KX,KY,KZ] = ndgrid(kk,kk,kk);
KR = sqrt(KX.^2 + KY.^2 + KZ.^2);

% Exclude DC and keep a thin spherical shell
keep_shell = (KR >= k_lo) & (KR <= k_hi) & (KR > 0);

% --- Optional anisotropy via cones about x,y,z ---------------------------
% A k-vector is kept if it lies within any chosen cone half-angle.
if any(coneDeg < 90)
    U = zeros([size(KX),3]);
    U(:,:,:,1) = KX ./ max(KR,eps);
    U(:,:,:,2) = KY ./ max(KR,eps);
    U(:,:,:,3) = KZ ./ max(KR,eps);

    coneKeep = false(size(KR));
    for a = 1:3
        if coneDeg(a) >= 90, continue; end
        coneKeep = coneKeep | (U(:,:,:,a) >= cosd(coneDeg(a)));
    end
    keep_shell = keep_shell & coneKeep;
end

% --- Sample nModes k-vectors uniformly from candidate shell --------------
idx = find(keep_shell);
if isempty(idx)
    error('No k-vectors passed shell+cone filters. Relax bandwidth/cones or change lambda_vox.');
end
if numel(idx) > nModes, idx = idx(randperm(numel(idx), nModes)); end

kxs = KX(idx); kys = KY(idx); kzs = KZ(idx);
phs = 2*pi*rand(size(kxs));                % random phases
amps = ones(size(kxs));                    % flat spectrum (OK for spinodoids)

% --- Build the periodic field (cosine sum = real, periodic) --------------
phi = zeros(N,N,N,'double');
s2p = 2*pi / L;
for m = 1:numel(kxs)
    phi = phi + amps(m) * cos(s2p*(kxs(m)*X + kys(m)*Y + kzs(m)*Z) + phs(m));
end

% Normalize to zero-mean, unit-std (helps thresholding)
phi = phi - mean(phi(:));
phi = phi ./ std(phi(:) + eps);

meta = struct('N',N,'L',L,'lambda_vox',lambda_vox,'k_target_idx',k_target_idx, ...
              'bandwidth',bandwidth,'nModes',numel(kxs),'coneDeg',coneDeg);
end
