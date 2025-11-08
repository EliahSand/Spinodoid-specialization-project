function [maskPad, rPad] = padRadial_stent(mask, rVec)
%PADRADIAL_STENT mirror one layer beyond inner/outer radii
dr = spacing_stent(rVec);
rPad = [max(rVec(1) - dr, 0), rVec(:).', rVec(end) + dr];
maskPad = false(size(mask) + [2 0 0]);
maskPad(2:end-1,:,:) = mask;
maskPad(1,:,:) = mask(1,:,:);
maskPad(end,:,:) = mask(end,:,:);
end

function [maskPad, zPad] = padAxial_stent(mask, zVec)
dz = spacing_stent(zVec);
zPad = [zVec(1) - dz, zVec(:).', zVec(end) + dz];
maskPad = false(size(mask) + [0 0 2]);
maskPad(:,:,2:end-1) = mask;
maskPad(:,:,1) = mask(:,:,1);
maskPad(:,:,end) = mask(:,:,end);
end

function [maskWrap, thetaWrap] = wrapTheta_stent(mask, thetaVec)
dtheta = spacing_stent(thetaVec);
thetaWrap = [thetaVec(:).', thetaVec(end) + dtheta];
maskWrap = cat(2, mask, mask(:,1,:));
end

function [X, Y, Z] = cylindricalGrid_stent(rVec, thetaVec, zVec)
[R, Theta, Z] = ndgrid(rVec, thetaVec, zVec);
X = R .* cos(Theta);
Y = R .* sin(Theta);
end

function d = spacing_stent(vec)
vec = vec(:).';
if numel(vec) > 1
    d = mean(diff(vec));
else
    d = max(vec, 1);
end
end
