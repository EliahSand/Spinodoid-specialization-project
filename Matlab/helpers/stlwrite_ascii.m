function stlwrite_ascii(filename, V, F)
fid = fopen(filename,'w');
if fid < 0, error('Cannot open %s', filename); end
fprintf(fid,'solid spinodoid\n');
Nf = faceNormals(V,F);
for i = 1:size(F,1)
    fprintf(fid,' facet normal %.7g %.7g %.7g\n', Nf(i,1), Nf(i,2), Nf(i,3));
    fprintf(fid,'  outer loop\n');
    for j = 1:3
        v = V(F(i,j),:);
        fprintf(fid,'   vertex %.7g %.7g %.7g\n', v(1), v(2), v(3));
    end
    fprintf(fid,'  endloop\n endfacet\n');
end
fprintf(fid,'endsolid spinodoid\n');
fclose(fid);
end
