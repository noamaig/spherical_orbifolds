%==========================================================================
% Example script for computing embeddings of sphere-topology meshes into
% sphere-topology spherical orbifolds.
%==========================================================================
%% init and load some data
init;
load  my_man_maxp.mat
%% embedding the mesh into the orbifold
% each orbifold type according to Fig. 2 in the paper, left 5
% examples (sphere-topology), choose a different orbifold_type to see
% embeddings into the different orbifolds.

orbifold_type=4;


switch orbifold_type
    case 1
        cones=[2 3 4];
    case 2
        cones=[2 3 3];
    case 3
        cones=[2 3 5];
    case 4
        k=4; % this can be any integer >1
        cones=[k 2 2];
    case 5
        inds=inds([1 3]); %this orbifold only has two cones so taking two selected vertices
        k=5; % this can be any integer >1
    otherwise
        error('orbifold type should be 1<=int<=5');
end
e=embed_from_data(V,T,cones,inds);

%% Drawing stuff
figure(1);
clf;

% draw the tiling of the embedding on the sphere
subplot(1,2,2);
e.draw('tileboundary',true,'colormap','bronze','tilepale',0.8);
ax=gca;
ax.CameraViewAngle=8.3;
campos([      -2.7762    8.5443    5.3982]);
camup([    0.1592   -0.4898    0.8572]);

% draw the source mesh with cuts
subplot(1,2,1);
e.draw('coords','source','colormap','bronze');
ax=gca;
ax.CameraViewAngle =  5.2845;
