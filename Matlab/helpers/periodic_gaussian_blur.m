function phi_s = periodic_gaussian_blur(phi, sigma_vox)
% Periodic Gaussian blur on a 3-D volume (not necessarily cubic) via FFT.
sz = size(phi);
assert(numel(sz)==3, 'periodic_gaussian_blur expects a 3-D volume.');
Nx = sz(1); Ny = sz(2); Nz = sz(3);
Phi = fftn(phi);
kkx = ifftshift(-floor(Nx/2):floor((Nx-1)/2));
kky = ifftshift(-floor(Ny/2):floor((Ny-1)/2));
kkz = ifftshift(-floor(Nz/2):floor((Nz-1)/2));
[KX,KY,KZ] = ndgrid(kkx,kky,kkz);
KR2 = (KX/Nx).^2 + (KY/Ny).^2 + (KZ/Nz).^2;
alpha = (2*pi*sigma_vox)^2;
G = exp(-alpha * KR2);
phi_s = real(ifftn(Phi .* G));
end
