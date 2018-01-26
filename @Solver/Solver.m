classdef Solver <handle
    properties
        inds;
        W;        
        I;
        J;
        adj;
        fix_gradient_tol;
        full_len_X;
        default_X;
        Wmat;
        precondition=true;
        radius=1;
        boundary;
    end
    methods
        function obj=Solver(adj,boundary,varargin)
            obj.boundary=boundary;
            parser = inputParser;
            parser.addOptional('fix_gradient_tol',0,@isnumeric);
            parser.addOptional('fix_inds',[],@(X)all(isnumeric(X)));
            parser.addOptional('default_x',zeros(length(adj),2));
            parser.parse(varargin{:});
            obj.full_len_X=length(adj);
            obj.default_X=parser.Results.default_x;
            if isempty(obj.default_X)
                obj.default_X=zeros(obj.full_len_X,2);
            end
            obj.Wmat=adj;
            obj.adj=adj~=0;
            [obj.I,obj.J]=find(adj);
            od=obj.I~=obj.J;
            obj.I=obj.I(od);
            obj.J=obj.J(od);
            W=full(adj(sub2ind(size(adj),obj.I,obj.J)));
            
            clamp=1e-3;
            negs=nnz(W<0);
            if negs
                %warning('clamping %d negative weights',negs);
            end
            W(W<0)=clamp;
            obj.W=double(single(W));
            
            Wmat_nonneg = sparse(obj.I,obj.J,obj.W,size(adj,1),size(adj,2));
            Wmat_nonneg = Wmat_nonneg - sparse(1:length(adj),1:length(adj),sum(Wmat_nonneg),size(adj,1),size(adj,2));
            obj.Wmat = Wmat_nonneg;
        end
    end
end

