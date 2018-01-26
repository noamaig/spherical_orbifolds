function x=computeFlattening(A,b,L)
%Euclidean tutte embedding - used for initialization

% L=V2A'*V2A;

% minimize ||Cx+d|| s.t. Ax=b
% minimize x^T L x s.t. Ax=b

n_vars = size(A,2);
n_eq = size(A,1);
M=[L A'; A sparse(n_eq,n_eq)];
rhs=[zeros(n_vars,1); b];
warning('off','MATLAB:nearlySingularMatrix')
x_lambda = M \ rhs;
warning('on','MATLAB:nearlySingularMatrix')
e=max(abs(M*x_lambda-rhs));
fprintf('error: %e\n',e);
if e>1e-1
    error('linear system not solved!');
end
if e>1e-6
    warning('linear system not solved!');
end

x = x_lambda(1:n_vars);
end

