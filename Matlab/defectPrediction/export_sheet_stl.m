function export_sheet_stl(caseDir, sheetMask, Svox, Lz)
%EXPORT_SHEET_STL Write the final defective sheet as an STL for visualization.
%
%   Vertices are scaled to metres (x,y by Svox; z by Lz/Nz).
%   Requires MATLAB R2018b+ (built-in stlwrite on triangulation objects).

stlPath = fullfile(caseDir, 'sheet.stl');

% sheetMask contains both solid base and spinodoid spin layers — full RVE.
fv = isosurface(double(sheetMask), 0.5);
if isempty(fv.faces)
    warning('export_sheet_stl:empty', 'No isosurface found for %s — STL not written.', caseDir);
    return;
end

Nz = size(sheetMask, 3);
dz = Lz / Nz;

% isosurface column layout: [row(y), col(x), z]
% Map back to physical coords and return as [x, y, z] for the STL.
verts = fv.vertices;
x_phys = (verts(:,2) - 1) * Svox;
y_phys = (verts(:,1) - 1) * Svox;
z_phys = (verts(:,3) - 1) * dz;
fv.vertices = [x_phys, y_phys, z_phys];

TR = triangulation(fv.faces, fv.vertices);
stlwrite(TR, stlPath);
end
