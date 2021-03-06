classdef OriginalBoundaryPiece < BoundaryPiece
    % a boundary segment of the disk orbifold
    
    properties
        %the reflection matrix mapping the boundary to itself
        A;
        %the normal to the supporting plane
        normal;
    end
    methods 
        function [obj]=OriginalBoundaryPiece(inds,spos,epos,color,normal)
            obj = obj@BoundaryPiece(inds,inds,spos,epos,color);
            if nargin>4
                obj.normal=normal;
            else
                obj.normal=cross(spos,epos);
                obj.normal=obj.normal/norm(obj.normal);
            end
            R=rotation_p2p(obj.normal,[1 0 0]);
            obj.A=R'*diag([-1 1 1]) *R;
        end
        function c=generateBoundaryConditions(obj)
             
            c=ReflectionBoundaryConditions(obj.normal,obj.inds(2:end-1));
        end
        function A=transformation(obj)
            A=obj.A;
        end
        function  N=supportingPlane(obj)
            N=obj.normal;
        end

    end
    
end

