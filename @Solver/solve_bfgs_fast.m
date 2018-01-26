
function [x,xinit,XinitFlat,log]=solve_bfgs_fast(obj)

lbfgs_memory = 3;

options =[];
options.MaxIter = 10000;
options.TolFun=0; %1e-10;
options.TolX=1e-10; 

% initialize
t_start = tic;
L = kron(obj.Wmat,speye(2));
[X,XinitFlat]=obj.boundary.initialSolution(L);
xinit=X;
X=X';X=X(:);
[A,b]=obj.boundary.generateBoundaryEquations(length(obj.adj));
log.t_init_flat = toc(t_start);

% setup lbfgs solver
t_start = tic;
fun = @(Y)obj.objective_and_grad(Y);
Ab = AffineSpace(A,b);
x0 = Ab.projectOnto(X);
% P = Precond_Identity();
P = Precond_Fixed(-Ab.N'*kron(obj.Wmat,eye(3))*Ab.N);
lbfgs_solver = OptimSolverLBFGS_NEW(fun, x0, Ab, P, lbfgs_memory);
log.t_init_lbfgs = toc(t_start);

% solve
t_start = tic;
%x = lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);
lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);
lbfgs_solver.updatePrecond(Precond_Identity());
lbfgs_solver.updateMemory(0);
x = lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);
log.t_lbfgs = toc(t_start);

% project back onto the sphere
fprintf('Vertex of smallest radius %g\n',min(sqrt(sum(reshape(x,3,[]).^2))));
fprintf('Vertex of largest radius %g\n',max(sqrt(sum(reshape(x,3,[]).^2))));
x=reshape(x,3,[]);
x=bsxfun(@mrdivide,x,sqrt(sum(x.^2)));
x = colStack(x);

% check boundary constraints
assert(max(abs((A*x-b)))<1e-6);
X=x(1:3:end);
            Y=x(2:3:end);
            Z=x(3:3:end);
            flatV=[X Y Z];
good=obj.boundary.validateSolution(flatV);
if max(good)>1e-6
    max(good)
    error();
end
end