classdef DiskOrbifoldStructure < EmbeddingStructure
    %implementation of the Orbifold structure for the disk-orbifold case
    
    properties
        V;
        T
        poses;
        inds;
        binds;
        cinds;
        ang;
    end
    methods
        function obj=DiskOrbifoldStructure(V,T,angs_or_poses,inds)
            tri=triangulation(T,V);
            obj.binds=tri.freeBoundary();
            obj.binds=obj.binds(:,1);
            assert(all(ismember(inds,obj.binds)),'all cone vertices must lie on the boundary.');
            firstInd=find(obj.binds==inds(1));
            obj.binds=[obj.binds(firstInd:end) ;obj.binds(1:firstInd-1)];
            [~,order]=ismember(inds,obj.binds);
            obj.cinds=order;
            assert(all(diff(order)>0),'cone vertices should be ordered by clockwise direction');
            if length(inds)==2
                assert(length(angs_or_poses)==1,'in case there are only two cones the 3rd argument should be the (single!) index of the cones.');
                assert(round(angs_or_poses)==angs_or_poses,'the index of the cone must be an interger.');
                poses=[0 0 1;
                    0 0 -1];
                obj.ang=2*pi/angs_or_poses;
                
            elseif length(inds)==3
                if length(angs_or_poses(:))==3
                    angs=angs_or_poses;
                    assert(length(angs)==length(inds),'number of cone vertices must match number of cone angles');
                    assert(all(round(angs_or_poses)==angs_or_poses),'the indices of the cone must be an interger');
                    poses=three_point_tile(angs,true);
                else
                    assert(all(size(angs_or_poses)==[3 3]),'if you give 3 cone inds and NOT 3 cone orders, the argument should be the positions of the 3 cones.');
                    poses=angs_or_poses;
                end
            else
                error('there can only be 2 or 3 cones!');
            end
            obj.V=V;
            obj.T=T;
            obj.inds=inds;
            obj.poses=poses;
        end
        function [boundary,cutmesh]=generate(obj)
            P={};
            if length(obj.inds)==3
                cols=linspecer(3);
                for i=1:3
                    if i<3
                        inds=obj.binds(obj.cinds(i):obj.cinds(i+1));
                        pos=obj.poses(i:i+1,:);
                        
                    else
                        inds=[obj.binds(obj.cinds(i):length(obj.binds))' obj.binds(1)];
                        pos=obj.poses([i 1],:);
                    end
                    P{i}=OriginalBoundaryPiece(inds,pos(1,:),pos(2,:),cols(i,:));
                end
            else
                cols=linspecer(2);
                inds=obj.binds(obj.cinds(1):obj.cinds(2));
                pos=obj.poses(1:2,:);
                P{1}=OriginalBoundaryPiece(inds,pos(1,:),pos(2,:),cols(1,:),[1 0 0]);
                
                inds=[obj.binds(obj.cinds(2):length(obj.binds))' obj.binds(1)];
                pos=obj.poses([2 1],:);
                P{2}=OriginalBoundaryPiece(inds,pos(1,:),pos(2,:),cols(2,:),[cos(-obj.ang) sin(-obj.ang) 0]); 
            end
            boundary=OrbifoldBoundary(P);
            cutmesh=[];
            cutmesh.V=obj.V;
            cutmesh.T=obj.T;
            cutmesh.oV=obj.V;
            cutmesh.oT=obj.T;
            cutmesh.new2old=1:length(obj.V);
            cutmesh.old2new=num2cell(1:length(obj.V));
            
        end
    end
    
    
    
    
end