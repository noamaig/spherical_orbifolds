function [o ] = orientation( Vmapped,T )

Vmapped=bsxfun(@mrdivide,Vmapped',sqrt(sum(Vmapped'.^2)))';

% %A=[a b c;
% %   d e f;
% %   g h i
% 
% %det(A)=aei-afh+bfg-bdi+cdh-ceg
% 
% a=Vmapped(T(:,1),1);
% b=Vmapped(T(:,1),2);
% c=Vmapped(T(:,1),3);
% d=Vmapped(T(:,2),1);
% e=Vmapped(T(:,2),2);
% f=Vmapped(T(:,2),3);
% g=Vmapped(T(:,3),1);
% h=Vmapped(T(:,3),2);
% i=Vmapped(T(:,3),3);
% o=a.*e.*i-a.*f.*h+b.*f.*g-b.*d.*i+c.*d.*h-c.*e.*g;

nFaces = size(T,1);
o = zeros(nFaces,1);
for ii = 1:nFaces
    o(ii) = det(Vmapped(T(ii,:),:));
end
end

