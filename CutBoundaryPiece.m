classdef CutBoundaryPiece < BoundaryPiece
    %represents a boundary piece of a mesh that's been cut open - i.e., it
    %has a twin and will have a rotation correlating them
    properties
        %pointer to the left twin of this segment
        leftTwin;
        %pointer to the right twin of this segment
        rightTwin;
        R;
        
    end
    
    methods
        function [obj]=CutBoundaryPiece(inds,original_inds,spos,epos,color,R)
            obj = obj@BoundaryPiece(inds,original_inds,spos,epos,color);
            if nargin>5
                obj.R=R;
            end
            
        end
        function N=supportingPlane(obj)
            p1=obj.startPos;
            p2=obj.endPos;
            if norm(p1+p2)>1e-8
             p1=obj.startPos;
                p2=obj.endPos;
                %first bring the two points to lie on x-y plane
                N=cross(p1,p2)';
                N=N/norm(N);
            else%antipodal 
                %vector orthogonal to p1,p2
                v=null([p1 ;ones(1,3)]); %there is an assumption here that null() is deterministic
                if obj.isLeftSide()
                    N=v;
                else
                    N=obj.twin().R*v;
                end
            end
        end
        function setR(obj)
            %set the rotation of this boundary piece according to the piece
            %and its twin
            if ~obj.isLeftSide() || ~isempty(obj.R)
                return;
            end
            p1=obj.startPos;
            p2=obj.endPos;
            p=[p1; p2];
            q1=obj.twin().startPos;
            q2=obj.twin().endPos;
            q=[q1; q2];
            if isempty(obj.R)
                assert(norm(p1+p2)>1e-8,'if the two cones are antipodal, you must supply a rotation matrix!');
                % want R s.t. P*R'=Q; or in other words min ||R*P'-Q'|| s.t. R rotation
                % solution to this is procustes: UEV=Q'*P
                % R=UV';
                [U,~,V]=closest_rotation(q'*p);
                obj.R=U*V';
                assert(norm(p*obj.R'-q,'fro')<1e-10);
                assert(det(obj.R)>0);
            end
        end
        function c=generateBoundaryConditions(obj)
            if ~obj.isLeftSide() %right side doesn't have boundary conditions - they were already defined on right 
                assert(isempty(obj.R),'there should''nt be a rotation set for the right-side of the cut');
                c=[];
                return;
            end
            %if left side we set boundary conditions
            obj.setR();
            %end points do not require this boundary conditions
            inds1=obj.inds(2:end-1);
            inds2=obj.twin().inds(2:end-1);
            c=RotationBoundaryConditions(obj.R,inds1,inds2);
        end
        function t=twin(obj)
            %return the twin segment for this segment
            if ~isempty(obj.leftTwin)
                t=obj.leftTwin;
            else
                t=obj.rightTwin;
            end
        end
        function o=isLeftSide(obj)
            %return true iff this segment is the "left" segment, false iff
            %right segment
            o=~isempty(obj.leftTwin);
        end
        function A=transformation(obj)
            A=obj.R;
        end
    end
    
end


