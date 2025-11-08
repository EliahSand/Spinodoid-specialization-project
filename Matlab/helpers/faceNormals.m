function Nf = faceNormals(V,F)
A = V(F(:,2),:) - V(F(:,1),:);
B = V(F(:,3),:) - V(F(:,1),:);
Nf = cross(A,B,2);
Nf = Nf ./ (sqrt(sum(Nf.^2,2)) + eps);
end
