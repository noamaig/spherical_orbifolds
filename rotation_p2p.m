function [ R ] = rotation_p2p( p1,p2 )
%R s.t. Rp1=p2

%R1*p1=[1 0 0; ....]
p1=fixvec(p1);
p2=fixvec(p2);
R1=[p1' ;null(p1')'];
if det(R1)<0
    R1(3,:)=-R1(3,:);
end
% assert(norm(R1*p1-[1 0 0]')<1e-10);

%R1*[1 0 0]=p2;
R2=[p2 null(p2')];
if det(R1)<0
    R1(:,3)=-R1(:,3);
end
R=R2*R1;
assert(norm(R*p1-p2)<1e-10);
assert(norm(R*R'-eye(3),'fro')<1e-10);



end


function x=fixvec(x)
if size(x,2)~=1
    x=x';
end
x=x/norm(x);

end