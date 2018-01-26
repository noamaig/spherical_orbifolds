classdef FixedDiskStructure < EmbeddingStructure
    % Embedding Structure for embedding into a fixed boundary - classical
    % Tutte into a convex polygon only on a sphere.
    
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
        function obj=FixedDiskStructure(V,T,angs_or_poses,inds)
            tri=triangulation(T,V);
            obj.binds=tri.freeBoundary();
            obj.binds=obj.binds(:,1);
            assert(all(ismember(inds,obj.binds)),'all cone vertices must lie on the boundary');
            firstInd=find(obj.binds==inds(1));
            obj.binds=[obj.binds(firstInd:end) ;obj.binds(1:firstInd-1)];
            [~,order]=ismember(inds,obj.binds);
            obj.cinds=order;
            assert(all(diff(order)>0),'cone vertices should be ordered by clockwise direction');
            
                
                    
                    assert(size(angs_or_poses,1)==length(inds),'number of cone vertices must match number of cone angles');
                    
                    
                
                    poses=angs_or_poses;
                
            
            obj.V=V;
            obj.T=T;
            obj.inds=inds;
            obj.poses=poses;
        end
        function [boundary,cutmesh]=generate(obj)
            P={};
            
                cols=linspecer(length(obj.cinds));
                for i=1:length(obj.cinds)
                    if i<length(obj.cinds)
                        inds=obj.binds(obj.cinds(i):obj.cinds(i+1));
                        pos=obj.poses(i:i+1,:);
                        
                    else
                        inds=[obj.binds(obj.cinds(i):length(obj.binds))' obj.binds(1)];
                        pos=obj.poses([i 1],:);
                    end
                    P{i}=FixedBoundaryPiece(inds,pos(1,:),pos(2,:),cols(i,:));
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