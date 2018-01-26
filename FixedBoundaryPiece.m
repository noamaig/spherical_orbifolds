classdef FixedBoundaryPiece < BoundaryPiece
    %boundary piece which should be held to place (fixed boundary conditions) 
    % - this is used when computing "classical" spherical tutte into a 
    % convex spherical polygon 
    properties
        normal;
        pos;
    end
    methods 
        function [obj]=FixedBoundaryPiece(inds,spos,epos,color)
            %inds - vertices along the piece
            %spos - position of the first vertex
            %epos - position of the last vertex
            %color - of the segment
            obj = obj@BoundaryPiece(inds,inds,spos,epos,color);
            obj.normal=cross(spos,epos);
                obj.normal=obj.normal/norm(obj.normal)';
            obj.pos=obj.positionsForInitialization();
                  
        end
        function c=generateBoundaryConditions(obj)
             
            c=FixedBoundaryConditions(obj.inds(2:end-1),obj.pos(2:end-1,:));
        end
        function A=transformation(obj)
            A=eye(3);
        end
        
            
            
        function  N=supportingPlane(obj)
            N=obj.normal;
        end

    end
    
end

