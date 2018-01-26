classdef ReflectionBoundaryConditions <BoundaryConditions
    % boundary conditions requiring vertices to satisfy v=Rv where R is a
    % reflection matrix. Essentialy this means that v lies on the inifinite
    % plane that R maps to itself, defined by a normal.
    
    properties
        
        normal;
        inds;
        
    end
    
    methods
        function obj=ReflectionBoundaryConditions(normal,inds)
            if size(normal,1)~=1
                normal=normal';
            end
            if size(inds,1)==1
                inds=inds';
            end
            
            obj.normal=normal/norm(normal);
            
            obj.inds=inds;
            
          
        end
        function [A,b]=constraints(obj,NV)
            I=repmat((1:length(obj.inds))',1,3)'; 
            J=[obj.inds*3-2 obj.inds*3-1 obj.inds*3]';
            V=repmat(obj.normal,length(obj.inds),1)';
            I=I(:);
            J=J(:);
            V=V(:);
            A=sparse(I,J,V,max(I),NV*3);
            b=zeros(size(A,1),1);
           
             
        end
        
        function good=validSolution(obj,V)
            
            good=max(abs(V(obj.inds,:)*obj.normal'));
        end
    end
    
end

