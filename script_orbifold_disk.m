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

orbifold_type=1;

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
e=embed_from_data(V,T,cones*2,inds);


%% Drawing the embedding

figure(1);
clf
% drawing the embedding on the sphere
subplot(1,2,2);
e.draw('tileboundary',true,'tilepale',0.8);
ax=gca;
ax.CameraViewAngle=8.3;
campos([   1.8559   -0.0047   10.3154]);
camup([   0 1 0]);

% drawing the original source mesh
subplot(1,2,1);
e.draw('coords','source');
ax=gca;
ax.CameraViewAngle =  5.2845;
camup([-0.6437   -0.3280   -0.6914]);
campos([  24.7040   59.0710  -42.88844]);
