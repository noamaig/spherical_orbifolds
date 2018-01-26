%==========================================================================
% Example script for computing embeddings of disk-topology meshes into
% disk-topology spherical orbifolds.
%==========================================================================
%% Init + load data
init;
load  fandisk.mat % NOTE THIS IS THE FANDISK WITH THE BOTTOM REMOVED SO IT'S A DISK!!!!!

%% embedding the mesh into the orbifold
% each orbifold type according to Fig. 2 in the paper, left 5
% examples (sphere-topology), choose a different orbifold_type to see
% embeddings into the different orbifolds.

%% setup
init;
load fandisk.mat
tri=triangulation(T,V);
b=tri.freeBoundary();
b=b(:,1);
inc=round(length(b)/7);
inds=b(1:inc:end);
theta=linspace(0,2*pi,length(inds)+1);
theta=theta(1:end-1)';
fac=1;
pos=[cos(theta)*fac sin(theta)*fac ones(length(theta),1)];
pos=bsxfun(@mrdivide,pos',sqrt(sum(pos'.^2)))';
%% computing
e=embed_from_data(V,T,pos,inds);


%% Drawing the embedding

figure(1);
clf
% drawing the embedding on the sphere
subplot(1,2,2);
e.draw('tileboundary',true,'tilepale',0.8);
[Xs,Ys,Zs]=sphere(500);
hold on
surf(Xs*0.99,Ys*0.99,Zs*0.99,'edgecolor','none','facecolor','black','facealpha',0.1)
ax=gca;
ax.CameraViewAngle=8.3;
campos([-2.4982   -2.8088    9.7835]);
camup([   1 0 0]);
light
% drawing the original source mesh
subplot(1,2,1);
e.draw('coords','source');
ax=gca;
ax.CameraViewAngle =  5.2845;
camup([-0.6437   -0.3280   -0.6914]);
campos([  24.7040   59.0710  -42.88844]);
