classdef FixedBoundaryConditions <BoundaryConditions
    % Object representing fixed point-constraint (i.e., "classic" Tutte) 
    % boundary conditions.
    
    properties        
        poses;
        inds;        
    end
    
    methods
        function obj=FixedBoundaryConditions(inds,poses)            
            if size(inds,1)==1
                inds=inds';
            end
            obj.inds=inds;
            if size(poses,2)~=3
                poses=poses';
            end
            obj.poses=poses;
            
        end
        function [A,b]=constraints(obj,NV)
            temp=(1:length(obj.inds))';
            I=[temp temp+length(obj.inds) temp+2*length(obj.inds)];
            J=[obj.inds*3-2 obj.inds*3-1 obj.inds*3];
            V=ones(size(I));
            I=I(:);
            J=J(:);
            V=V(:);
            A=sparse(I,J,V,max(I),NV*3);
            b=obj.poses(:);              
        end        
        function good=validSolution(obj,V)
            
            good=max(max(abs(V(obj.inds,:)-obj.poses)));
        end
    end
    
end

