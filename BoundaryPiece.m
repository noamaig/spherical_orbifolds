classdef BoundaryPiece < handle
    %object representing "part" of the boundary, i.e. the line between two
    %cones.
    properties
        %pointer to prev piece
        prev;
        %next piece
        next;
        %indices of the vertices on this segment on the disk mesh (i.e., cut)
        inds;
        %indices of the vertices on this segment on the original mesh (i.e., uncut)
        originalInds;
        %pos of start point
        startPos;
        %pos of end point
        endPos;
        %color of the curve
        color;
    end
    methods
        function obj=BoundaryPiece(inds,originalInds,spos,epos,color)
            assert(~isempty(inds));
            obj.originalInds=originalInds;
            obj.inds=inds;
            obj.startPos=spos;
            obj.endPos=epos;
            obj.color=color;
        end
        function [s,e]=endPoints(obj)
            s=obj.inds(1);
            e=obj.inds(end);
        end
        function draw(obj,V)
            v=V(obj.inds,:);
            line(v(:,1),v(:,2),v(:,3),'linewidth',5,'color',obj.color);
        end
        function p=positionsForInitialization(obj)
            %positions of the vertices along the segment according to the
            %basic polygon (used for initialization with tutte to convex
            %polygon)
            p1=obj.startPos;
            p2=obj.endPos;
            %do slurp between the two end points
            N=obj.supportingPlane();
            if size(N,2)~=1
                N=N';
            end
            assert(all(size(N)==[3 1]));
            R1=rotation_p2p(N,[0 0 1]);
            assert(norm(R1*N-[0;0;1])<1e-8);
            %points rotated to xy plane
            q1=R1*p1';
            q2=R1*p2';
            q1=q1(1:2);
            q2=q2(1:2);
            
            %now rotate points so q1 is at (0,1)
            %q1=cos(t) sin(t)
            %we want R2*q1=[0 1];
            %R2=[-i*q1;q1;];
            
            %                 R2=[q1(2) -q1(1);q1'];
            R2=[q1';
                -q1(2) q1(1)];
            assert(norm(R2*q1-[1;0])<1e-8);
            a=R2*q2;
            %a is q2 and all we need to do is take its angle
            %                 assert(a(2)>=0);
            %q1 is at (1,0) so its angle on circle is pi/2
            theta1=0;
            %compute angle of q2
            theta2=atan2(a(2),a(1));
            if theta2>pi
                theta2=pi-theta2;
            end
            %interpolate angles to all vertices
            theta=linspace(theta1,theta2,length(obj.inds));
            %points from angles
            p=[cos(theta)' sin(theta)'];
            %reverse the transformations to bring them back to place
            %move q1 from (0,1) to original position
            p=p*R2;
            %back into 3d from 2d
            p=[p zeros(length(theta),1)];
            %move xy plane back to original plane
            p=p*R1;
            %check that both end points are at their original position
            assert(norm(p1-p(1,:))<1e-8);
            assert(norm(p2-p(end,:))<1e-8);
            
            
            
        end
    end
    
    methods (Abstract)
        %method for generarting boundary conditions for the vertices on
        %this segment
        c=generateBoundaryConditions(obj);
        %the automorphism associated with this segement
        %(rotation\reflection)
        A=transformation(obj);
        %normal of the plane supporting this piece
        N=supportingPlane(obj);
    end
    
end

