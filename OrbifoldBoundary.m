classdef OrbifoldBoundary < handle
    %object that holds all the boundary pieces
    
    properties
        boundaryPieces;
        constraints={};
        coneInds;
        conePos;
        origConeInds;
        coneCols;
    end
    
    methods
        function obj=OrbifoldBoundary(boundaryPieces)
            
            obj.boundaryPieces=boundaryPieces;
            obj.conePos=[];
            obj.coneInds=[];
            origConeInds=[];
            for i=1:length(boundaryPieces)
                p=boundaryPieces{i};
                [s,e]=p.endPoints();
                obj.coneInds=[obj.coneInds;s];
                obj.coneInds=[obj.coneInds;e];
                origConeInds=[origConeInds;reshape(p.originalInds([1 end]),2,1)];
                obj.conePos=[obj.conePos;boundaryPieces{i}.startPos;boundaryPieces{i}.endPos];
            end
            [I,J]=unique(obj.coneInds);
            obj.coneInds=I';
            obj.origConeInds=origConeInds(J);
            obj.conePos=obj.conePos(J,:);
            
            [I,~,J]=unique(obj.origConeInds);
            %             cols=linspecer(length(I));
            if length(I)<=3
                cols=[1 0.8 0;
                    0.4 1 0.7;
                    1 0 0.5;
                    
                    ];
            else
                cols=linspecer(length(I));
            end
            obj.coneCols=cols(J,:);
            for i=1:length(obj.boundaryPieces)
                p=obj.boundaryPieces{i};
                c=p.generateBoundaryConditions();
                if isempty(c)
                    continue;
                end
                obj.constraints{end+1}=c;
            end
            
            for i=1:length(obj.coneInds)
                
                s=obj.coneInds(i);
                pos=obj.conePos(i,:);
                obj.constraints{end+1}=FixedBoundaryConditions(s,pos);
            end
            
        end
        function [A,b]=generateBoundaryEquations(obj,NV)
            %generate the entire equation system for the boundary vertices
            % NV - integer, number of vertices in the mesh (needed so the
            % width of the constraint matrix A is correct).
            if isempty(obj.constraints)
                error('need to call setConstraints');
            end
            As=cell([1 length(obj.constraints)]);
            bs=cell([1 length(obj.constraints)]);
            for i=1:length(obj.constraints)
                [As{i},bs{i}]=obj.constraints{i}.constraints(NV);
            end
            A=cell2mat(As');
            b=cell2mat(bs');
        end
        function good=validateSolution(obj,V)
            %return a vector where good(i)==true iff i'th constraint is
            %satisfied
            if isempty(obj.constraints)
                error('need to call setConstraints');
            end
            good=zeros(length(obj.constraints),1);
            for i=1:length(obj.constraints)
                good(i)=max(obj.constraints{i}.validSolution(V));
            end
        end
        function [v0,flatV]=initialSolution(obj,L)
            %create an initial embedding, a starting point for the optimization
            % input: L - a laplacian matrix to use for embedding
            % output: v0 - the embedding on the sphere, flatV - the planar
            % embedding that was used (that was projected back on the
            % sphere)
            cons=PosConstraints(length(L)/2);
            p=obj.conePos;
            if size(p,1)==2
                init_p1=obj.boundaryPieces{1}.positionsForInitialization();
                init_p2=obj.boundaryPieces{2}.positionsForInitialization();
                
                p=[init_p1;init_p2];
                
            end
            
            N=findHemisphere(p);
            R=rotation_p2p(N,[0 0 1]);
            
            for i=1:length(obj.boundaryPieces)
                p=obj.boundaryPieces{i};
                inds=p.inds;
                init_p=p.positionsForInitialization();
                pp=init_p*R';
                assert(all(pp(:,3)>=0-eps));
                for j=1:length(inds)
                    
                    cons.addConstraint(inds(j),1,orthographic_to_plane(pp(j,:)));
                end
                
            end
            
            
            
            
            x=computeFlattening(cons.A,cons.b,L);
            
            X=x(1:2:end);
            Y=x(2:2:end);
            flatV=[X Y ];
            
            v0=orthographic_to_sphere(flatV);
            v0=v0*R;
        end
        function draw(obj,V)
            %draw the boundary
            for i=1:length(obj.boundaryPieces)
                p=obj.boundaryPieces{i};
                p.draw(V);
            end
            hold on;
            
            v=V(obj.coneInds,:);
            scatter3(v(:,1),v(:,2),v(:,3),150,'filled','CData',obj.coneCols);
        end
        function As=tilingTransformations(obj)
            %automorphisms of the orbifold
            t=Tiler();
            for i=1:length(obj.boundaryPieces)
                p=obj.boundaryPieces{i};
                A=p.transformation();
                if isempty(A)
                    continue;
                end
                t.addIfNew(A);
                t.addIfNew(A');
            end
            
            t.tile();
            As=t.As;
        end
        
    end
end
