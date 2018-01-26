function [os, isdisk ] = createOrbifoldStructure(V,T,angs_or_poses,inds)
%hi level func for generating the orbifold strucutre. Gets the original
%mesh (V,T) the angles or symmetry order of the cones (angs_or_poses) and
%the indices of the cone vertices (inds) and returns the correct orbifold
%structure (os) and a flag of whether the mesh is a disk or a sphere.
euler=euler_char(V,T);
if euler==1
    isdisk=true;
    if min(size(angs_or_poses))==1
        os=DiskOrbifoldStructure(V,T,angs_or_poses,inds);
    else
        os=FixedDiskStructure(V,T,angs_or_poses,inds);
    end
elseif euler==2
    isdisk=false;
    os=SphereOrbifoldStructure(V,T,angs_or_poses,inds);
else
    error('euler characteristic is wrong! %d',euler);
end

end
function euler=euler_char(V,T)
n = size(T,2);

e = nchoosek(1:n,2);
A = sparse(T(:,e(:,1)),T(:,e(:,2)),1,max(T(:)),max(T(:)));
[EI,EJ] = find(tril(A+A'));
E = [EJ EI];
ne=length(E);
nt=length(T);
nv=length(V);
euler=nv-ne+nt;
end