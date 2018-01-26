function [ Y ] = stereographic_uv( V )
V=normr(V);
[ n ] = findHemisphere( V );
R1=rotation_p2p(n,[0 0 1]);
R2=[1 0 0;
    0 0 -1;
    0 -1 0];
V=V*R1'*R2';
Y=[V(:,1)./(1-V(:,3)) V(:,2)./(1-V(:,3))];

end

