function P=three_point_tile(cone_order,disk)
%given 3 angles of the cones, return the 4 vertices of the tile

% angles 1 and 3 are divided because we take half a triangle, angle 2 is
% divided for another reason - that cone is a hybrid
if nargin<2
    disk=false;
end

angs=2*pi./cone_order;
if ~disk
    realangs=angs/2;
else
    realangs=angs;
end

%  1                 ^
%                    |
%                   z axis
%  3    2
%                    x axis ->
%
[a1, b1, c1, a2, b2, c2] = aaa(realangs(1), realangs(2), realangs(3));
t=[a1 b1 c1;a2 b2 c2];
for i=1:3
    if i==3
        error('no valid solution - probably your choice of cones doesn''t match any of the orbifolds?');
    end
    if ~any(isnan(t(i,:)))
        t=t(i,:);
        break;
    end
    
end

p1=[0 0 1]';
%p2=p1 rotated on y-z plane by angle_2
ang2=t(2);
%  R=[1 0 0;
%     0  cos(ang) sin(ang);
%     0  -sin(ang) cos(ang)]
p2=[0 sin(ang2) cos(ang2)]';

%first generate p3 on the y-z plane
ang3=t(3);
ptemp=[0 sin(ang3) cos(ang3)]';
%now rotate it according to ang1
ang=-realangs(1);
R=[cos(ang) sin(ang) 0;
    -sin(ang) cos(ang) 0;
    0 0 1;];
p3=R*ptemp;
if disk
    P=[p1';p2';p3'];
    return;
end
p4=R'*ptemp;
P=[p1';   
p3';
p2';
p4';];
end