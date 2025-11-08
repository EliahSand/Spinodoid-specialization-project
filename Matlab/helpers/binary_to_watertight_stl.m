function meshStats = binary_to_watertight_stl(solidMask, L, outPath)
% Create a closed surface of the solid region by marching cubes on the binary.
% - solidMask: logical NxNxN (true = solid)
% - L: physical box size

fv = isosurface(double(solidMask), 0.5);
if isempty(fv.vertices)
    error('Empty surface: the solid mask is all false or extremely sparse.');
end

N = size(solidMask,1);
S = L / N;
V = (fv.vertices - 1) * S;
F = double(fv.faces);

stlwrite_ascii(outPath, V, F);
meshStats = struct('numFaces', size(F,1), 'numVertices', size(V,1));
end
