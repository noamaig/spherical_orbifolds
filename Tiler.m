classdef Tiler<handle
    %object for using the automorphisms of the orbifold to generate a
    %tiling of the sphere
    
    
    properties
        %stack of transformations
        As=zeros(3,3,0);
    end
    methods
        function obj=Tiler()
            %first transformation is identiry
            obj.As(:,:,1)=eye(3);
        end
        function addIfNew(obj,addA)
            %add matrix if it's a new one
            if obj.isNew(addA)
                obj.As(:,:,end+1)=addA;
%                 fprintf('added trans to the tiling group, its size: %d\n',size(obj.As,3));
            end
            
        end
        function tile(obj)
            %generate the tiling
            i=1;
            while(i<=size(obj.As,3))
                %get the current transformation
                A=obj.As(:,:,i);
                %go over all other transformations
                for j=1:size(obj.As,3)
                    %compose the current trans with the prev. one to get a
                    %new
                    addA=obj.As(:,:,j)*A;
                    %add it if it's a new one
                    obj.addIfNew(addA);
                    %check from the other side also.
                    addA=A*obj.As(:,:,j);
                    obj.addIfNew(addA);
                end
                i=i+1;
            end
        end
        function new=isNew(obj,B)
            for j=1:size(obj.As,3)
                A=obj.As(:,:,j);
                if norm(A-B,'fro')<1e-8
                    new=false;
                    return
                end
            end
            new=true;
        end
    end
end

