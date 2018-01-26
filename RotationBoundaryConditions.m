classdef RotationBoundaryConditions <BoundaryConditions
    % boundary conditions requiring that two sets of vertices, v and w,
    % satisfy v=Rw where R is a rotation matrix.
    
    properties
        
        R;%the rotation
        inds1;%first set of vertex indices
        inds2;%second set of vertex indices
    end
    
    methods
        function obj=RotationBoundaryConditions(R,inds1,inds2)
            
            obj.R=R;
            if size(inds1,1)==1
                inds1=inds1';
            end
            if size(inds2,1)==1
                inds2=inds2';
            end
            obj.inds1=inds1;
            obj.inds2=inds2;
        end
        
        function [A,b]=constraints(obj,NV)
            %R*v...
            %R_11*v_x+R_12*v_y+R_13*v_z=w_x
            %3 rows of constraitns for each vertex each row with 4 vars
            
            %work with N-on-4 matrices and only in the end column-stack them,
            % to maintain sanity. 
            
            % [1 1 1 1...
            %  2 2 2 2
            %  . . . .
            % 3|V| 3|V| 3|V| 3|V|]
            I=repmat((1:(3*length(obj.inds1)))',1,4)'; 
            
            %[vx vy vz wx; %%(times |V|)
            % vx vy vz wy; %%(times |V|)
            % vx vy vz wz] %%(times |V|)
            leftSide=repmat([obj.inds1*3-2 obj.inds1*3-1 obj.inds1*3],3,1);
            rightSide=[obj.inds2*3-2;obj.inds2*3-1;obj.inds2*3];
            J=[leftSide rightSide]';
            
            %[ R(1,:) -1; %%(times |V|)
            %[ R(2,:) -1; %%(times |V|)
            %[ R(3,:) -1] %%(times |V|)
            V=[];
            for i=1:3
                line=[repmat(obj.R(i,:),length(obj.inds1),1) -ones(length(obj.inds1),1)];
                V=[V; line];
            end
            I=I(:);
            J=J(:);
            V=V';
            V=V(:);
            A=sparse(I,J,V,max(I),round(NV)*3);
            b=zeros(size(A,1),1);
           
             
        end
        function good=validSolution(obj,V)
            
            good=sum((V(obj.inds1,:)*obj.R'-V(obj.inds2,:)).^2,2);
        end
    end
    
end

