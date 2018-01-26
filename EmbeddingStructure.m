classdef EmbeddingStructure < handle
    %the object that generates the entire structure of the embedding, gets
    %the mesh & the cones, and outputs the disk mesh of the basic tile and
    %the boundary object which is responsible for generating
    %boundary conditions etc.
    
    properties
    end
    
    
    methods (Abstract)
        
        [boundary,cutmesh]=generate(obj);
    end
end

