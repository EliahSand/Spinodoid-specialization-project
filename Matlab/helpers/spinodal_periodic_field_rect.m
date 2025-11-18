function [phi, meta] = spinodal_periodic_field_rect(Nvec, Lvec, lambda_vox, bandwidth, nModes, coneDeg)
% Rectangular periodic Gaussian random field via discrete Fourier modes.
% Nvec: [Nx Ny Nz] grid samples along each axis
% Lvec: [Lx Ly Lz] physical lengths along each axis

Nvec = round(Nvec(:)');
assert(numel(Nvec)==3, 'Nvec must have three elements [Nx Ny Nz].');
Lvec = Lvec(:)';
assert(numel(Lvec)==3, 'Lvec must have three elements [Lx Ly Lz].');

Nx = Nvec(1); Ny = Nvec(2); Nz = Nvec(3);
Lx = Lvec(1); Ly = Lvec(2); Lz = Lvec(3);

dx = Lx / Nx;
dy = Ly / Ny;
dz = Lz / Nz;
vox_mean = mean([dx dy dz]);
lambda_phys = lambda_vox * vox_mean;
k_target = 2*pi / lambda_phys;
k_lo = (1 - bandwidth) * k_target;
k_hi = (1 + bandwidth) * k_target;

kkx = -floor(Nx/2):floor((Nx-1)/2);
kky = -floor(Ny/2):floor((Ny-1)/2);
kkz = -floor(Nz/2):floor((Nz-1)/2);
[KX,KY,KZ] = ndgrid(kkx,kky,kkz);

Kphys = 2*pi * sqrt( (KX./max(Lx,eps)).^2 + (KY./max(Ly,eps)).^2 + (KZ./max(Lz,eps)).^2 );
keep_shell = (Kphys >= k_lo) & (Kphys <= k_hi) & ((KX~=0) | (KY~=0) | (KZ~=0));

if any(coneDeg < 90)
    normK = sqrt((KX./max(Lx,eps)).^2 + (KY./max(Ly,eps)).^2 + (KZ./max(Lz,eps)).^2);
    Ux = (KX./max(Lx,eps)) ./ max(normK, eps);
    Uy = (KY./max(Ly,eps)) ./ max(normK, eps);
    Uz = (KZ./max(Lz,eps)) ./ max(normK, eps);
    cone = false(size(normK));
    if coneDeg(1) < 90, cone = cone | (Ux >= cosd(coneDeg(1))); end
    if coneDeg(2) < 90, cone = cone | (Uy >= cosd(coneDeg(2))); end
    if coneDeg(3) < 90, cone = cone | (Uz >= cosd(coneDeg(3))); end
    keep_shell = keep_shell & cone;
end

idx = find(keep_shell);
if isempty(idx)
    error('spinodal_periodic_field_rect:noModes', ...
        'No Fourier modes satisfy the shell/cone filters. Adjust bandwidth or lambda_vox.');
end
if numel(idx) > nModes
    idx = idx(randperm(numel(idx), nModes));
end

kxs = KX(idx);
kys = KY(idx);
kzs = KZ(idx);
phs = 2*pi*rand(size(kxs));
amps = ones(size(kxs));

x = (0:Nx-1)/Nx * Lx;
y = (0:Ny-1)/Ny * Ly;
z = (0:Nz-1)/Nz * Lz;
[X,Y,Z] = ndgrid(x,y,z);

phi = zeros(Nx,Ny,Nz,'double');
for m = 1:numel(kxs)
    arg = 2*pi * (kxs(m)*X/Lx + kys(m)*Y/Ly + kzs(m)*Z/Lz) + phs(m);
    phi = phi + amps(m) * cos(arg);
end

phi = phi - mean(phi(:));
phi = phi ./ (std(phi(:)) + eps);

meta = struct('Nvec',Nvec,'Lvec',Lvec,'lambda_vox',lambda_vox, ...
              'bandwidth',bandwidth,'nModes',numel(kxs),'coneDeg',coneDeg, ...
              'k_target',k_target);
end
