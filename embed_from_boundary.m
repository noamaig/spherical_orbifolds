function [ e ] = embed_from_boundary( boundary,cut_mesh )
%hi level function that computes the embedding from the given boundary
%object and the cut-open mesh
tid=tic;
V=cut_mesh.V;
T=cut_mesh.T;
L=cotmatrix(V,T);
s=Solver(L,boundary);
timing.pre=toc(tid);
fprintf('*** Pre processing: ');
toc(tid);
tid=tic;
[x,xinit,XinitFlat,timing.lbfgs]=s.solve_bfgs_fast();
timing.solve=toc(tid);
fprintf('*** solve: ');
toc(tid);
tid=tic;
X=[x(1:3:end) x(2:3:end) x(3:3:end)];
e=Embedding(cut_mesh,X,boundary,xinit,XinitFlat);
timing.post=toc(tid);
fprintf('*** Post processing: ');
toc(tid);
e.timing=timing;
end

