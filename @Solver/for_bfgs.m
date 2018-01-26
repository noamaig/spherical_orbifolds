function [o,g]=for_bfgs(obj,X)

% o=0;
% g=zeros(size(X));
% return
%             if any(d>1)
%                 o=inf;
%                 g=zeros(length(X)*2,1);
%                 return;
%             end

%             obj.check_derivs(X);
%             o=obj.objectives(X);
%             j=obj.GN_jacobian(X);
%             g=j'*o;

[o,g]=obj.objective_and_grad(X);

Y=reshape(X,[3 length(X)/3 ])';

%f(x)=(x^2+y^2+z^2-1)^2
%df=4(x^2^y^2+z^2-1)x
sigma=10;%1000;
f=sum(Y.^2,2)-1;

Obar=sigma*sum(f.^2);
Gbar=sigma*4*repmat(f,1,3).*Y;

% Obar=sigma*exp(sum(f.^2)-1);
% Gbar=sigma*4*exp(sum(f.^2))*repmat(f,1,3).*Y;
Gbar=Gbar';
Gbar=Gbar(:);
o=o+Obar;
g=g+Gbar;
%next step is heuristic division by the metric for
%"precondtioning"

% disp(Obar)



end
