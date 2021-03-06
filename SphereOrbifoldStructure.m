classdef SphereOrbifoldStructure < EmbeddingStructure
    % the embedding structure for a spherical orbifold.
    
    properties
        V;
        T
        poses;
        inds;
        R;
        
    end
    
    methods
        function obj=SphereOrbifoldStructure(V,T,angs_or_poses,inds)
            %V,T - the uncut spherical mesh
            %angs_or_poses - the cones' symmetry order or the angles of the cones
            %inds - the indices of the vertices that will be cones
            if length(inds)==2
                assert(length(angs_or_poses)==1,'in case there are only two cones the 3rd argument should be the (single!) index of the cones');
                assert(round(angs_or_poses)==angs_or_poses,'the index of the cone must be an interger');
                poses=[0 0 1;
                    0 0 -1];
                ang=2*pi/angs_or_poses;
                obj.R=[cos(ang) sin(ang) 0;
                    -sin(ang) cos(ang) 0;
                    0 0 1];
            elseif length(inds)==3
                if length(angs_or_poses(:))==3
                    angs=angs_or_poses;
                    assert(length(angs)==length(inds),'number of cone vertices must match number of cone angles');
                    assert(all(round(angs_or_poses)==angs_or_poses),'the indices of the cone must be an interger');
                    poses=three_point_tile(angs);
                else
                    poses=angs_or_poses;
                    assert(size(poses,1)==length(inds)*2-2,'number of cone angles must equal number of cones');
                end
            else
                error('there can only be 2 or 3 cones!');
            end
            obj.V=V;
            obj.T=T;
            obj.inds=inds;
            obj.poses=poses;
            
        end
        function [boundary,cutter]=generate(obj)
            fixedPairs=[1:length(obj.inds)-1;2:length(obj.inds)]';
            tree=sparse(fixedPairs(:,1),fixedPairs(:,2),1,length(obj.inds),length(obj.inds));
            cutter=TreeCutter(obj.V,obj.T,tree,obj.inds,1);%+(rand(size(obj.M_orig.V'))*20)*diag([1 100 100])+rand(size(obj.M_orig.V'))*20
            cutter.cutTree();
            if size(obj.poses,1)>2
                leftSide=2:(length(obj.poses)/2);
                rightSide=(length(obj.poses)):-1:(length(obj.poses)/2+2);
                leftposes=obj.poses([1 leftSide length(obj.poses)/2+1],:);
                rightposes=obj.poses([1 rightSide length(obj.poses)/2+1],:);
            else
                leftposes=obj.poses;
                rightposes=obj.poses;
            end
            P={};
            Q={};
            cols=linspecer(length(cutter.pathPairs));
            for i=1:length(cutter.pathPairs)
                
                P{end+1}=CutBoundaryPiece(cutter.pathPairs{i}(:,1),cutter.new2old(cutter.pathPairs{i}(:,1)),leftposes(i,:),leftposes(i+1,:),cols(i,:),obj.R);
                Q{end+1}=CutBoundaryPiece(cutter.pathPairs{i}(:,2),cutter.new2old(cutter.pathPairs{i}(:,1)),rightposes(i,:),rightposes(i+1,:),cols(i,:));
                P{end}.leftTwin=Q{end};
                Q{end}.rightTwin=P{end};
            end
            
            P=horzcat(P,Q);
            boundary=OrbifoldBoundary(P);
            
        end
    end
    
end

