classdef BoundaryConditions < handle
    %abstract class for the boundary conditions object
    
    properties
    end
   
    methods (Abstract)
        
        good=validSolution(V); %return true iff solution is valid 
        [A,b]=constraints(obj,NV); % generate the constraint matrix and RHS
    end
        
   
    
end

