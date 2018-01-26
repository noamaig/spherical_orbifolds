classdef Embedding < handle
    %Object representing the computed embedding of the mesh into the orbifold
    
    properties
        
        
        Y;
        cut_mesh;
        boundary;
        smin;
        smax;
        V2A;
        areas;
        As;
        orientation;
        flipped;
        frobenius;
        spectralEmbedding;
        timing;
        xinit;
        xinitflat;
    end
    
    methods
        function e=rerun(obj)
            %recomputre the embedding, return it without changing this
            %embedding
            e= embed_from_boundary(obj.boundary,obj.cut_mesh);
        end
        function v=T(obj)
            %triangles of the basic tile
            v=obj.cut_mesh.T;
        end
        
        
        function obj=Embedding(disk_mesh,Y,boundary,xinit,xinitflat)
            % disk_mesh - the disk topolgy mesh for the basic tile, either
            % a cut-open mesh, or a mesh that was originally a disk.
            % Y - the positions of the vertices
            % boundary - the OrbifoldBoundary object
            % xinit - the initialization of the embedding before
            % optimization on the sphere.
            % xinitflat - the Eculidean Tutte embedding that was projected
            % onto the sphere as the initialization
            obj.xinitflat=xinitflat;
            obj.xinit=xinit;
            obj.cut_mesh=disk_mesh;
            obj.Y=Y;
            obj.boundary=boundary;
            
            L=cotmatrix(obj.cut_mesh.oV,obj.cut_mesh.oT);
            [E,S]=eigs(L,4,'sm');
            S=diag(S);
            [~,ind]=min(abs(S));
            E(:,ind)=[];
            obj.spectralEmbedding=E(obj.cut_mesh.new2old,:);
            obj.orientation=orientation(obj.Y,obj.cut_mesh.T);
            
            obj.flipped=obj.orientation<0;
            if(any(obj.flipped))
                disp('*********************************');
                fprintf('*** THERE ARE %d FLIPPED TRIS!!! ***\n',nnz(obj.flipped));
                disp('*********************************');
            end
        end
        
        
        function  meshcol=draw(obj,varargin)
            % draw the mesh\embedding, see the output() function for
            % arguments
            [~,meshcol]=obj.output(varargin{:});
        end
        
        function [S,meshcol]=output(obj,varargin)
            %varargin - see below
            p = inputParser;
            p.addOptional('edges',false,@islogical);
            p.addOptional('edgealpha',1,@(x)(isnumeric(x)&&x>=0&&x<=1));
            p.addOptional('coords','target');
            p.addOptional('tile',true,@islogical);
            p.addOptional('boundary',true,@islogical);
            p.addOptional('color','vcol');
            p.addOptional('colorbound',5,@(x)(isnumeric(x)&&x>=1));
            p.addOptional('colormap','hsv',@isstr);
            p.addOptional('colormodel','eig',@isstr);
            p.addOptional('tileboundary',false,@islogical);
            p.addOptional('tilepale',0,@(x)(isnumeric(x)&&x>=0&&x<=1));
            p.addOptional('normalshading',true,@islogical);
            p.addOptional('normalspecular',false,@islogical);
            p.addOptional('drawflips',true,@islogical);
            p.addOptional('colorchannels',[]);
            p.addOptional('uv',false,@islogical);
            p.parse(varargin{:});
            onSphere=true;
            if ischar(p.Results.coords)
                if strcmp(p.Results.coords,'source')
                    X=obj.cut_mesh.V;
                    onSphere=false;
                elseif strcmp(p.Results.coords,'target')
                    X=obj.Y;
                elseif strcmp(p.Results.coords,'init')
                    X=obj.xinit;
                elseif strcmp(p.Results.coords,'flatinit')
                    X=obj.xinitflat;
                    
                    onSphere=false;
                else
                    error('unknown coords: %s',p.Results.coords);
                end
            else
                X=p.Results.coords;
            end
            if size(X,2)==2
                X=[X zeros(length(X),1)];
            end
            if p.Results.tile && onSphere
                tileAs=obj.boundary.tilingTransformations();
            else
                tileAs=[];
            end
            if all(X(:,3)==0)
                
                BBmin=min(X);
                BBmax=max(X);
                X=bsxfun(@minus,X,mean(X));
                X=2*X/max(BBmax-BBmin);
            end
            if strcmp(p.Results.colormodel,'eig')
                vcol=obj.spectralEmbedding;
                for i=1:size(vcol,2)
                    if vcol(1,i)<0
                        vcol(:,i)=-vcol(:,i);
                    end
                end
                if size(vcol,2)==2
                    vcol=[vcol zeros(length(vcol),1)];
                end
            elseif strcmp(p.Results.colormodel,'rgb')
                vcol=obj.cut_mesh.V;
            else
                tri=triangulation(obj.cut_mesh.oT,obj.cut_mesh.oV);
                vcol=(tri.vertexNormal());
                vcol=vcol(obj.cut_mesh.new2old,:);
            end
            
            if ~(strcmp(p.Results.colormodel,'eig')||strcmp(p.Results.colormodel,'rgb'))
                vcol=abs(vcol);
                for i=1:3
                    vcol(:,i)=vcol(:,i)-min(vcol(:,i));
                    vcol(:,i)=vcol(:,i)./max(vcol(:,i));
                end
                vcol(isnan(vcol))=0;
                vcol=bsxfun(@rdivide,vcol,sum(vcol,2));
                col=zeros(size(vcol));
                if strcmp(p.Results.colormodel,'gold')
                    ccol=[1 0.9 0;0.9 0.9 0.9;[139,69,19]/255];
                    
                else
                    
                    ccol=p.Results.colorchannels;
                    if isempty(ccol)
                        ccol=ones(3)-eye(3);
                    end
                end
                for i=1:3
                    
                    col=col+vcol(:,i)*ccol(i,:);
                end
                vcol=col;
            else
                for i=1:3
                    vcol(:,i)=vcol(:,i)-min(vcol(:,i));
                    vcol(:,i)=vcol(:,i)./max(vcol(:,i));
                end
                vcol(isnan(vcol))=0;
            end
            
            
            
            
            
            tilecol=vcol*(1-p.Results.tilepale)+p.Results.tilepale;
            V=obj.cut_mesh.V();
            if any(V(:,3)~=0)
                tri=triangulation(obj.cut_mesh.oT,obj.cut_mesh.oV);
                normals=(tri.vertexNormal());
                normals=normals(obj.cut_mesh.new2old,:);
                normals=abs(normals*[1 0 0]');
                normals=normals-min(normals);
                normals=normals/max(normals);
                t=0.4;
                normals=normals*t+(1-t);
                if p.Results.normalshading
                    vcol=vcol.*repmat(normals,1,3);
                    tilecol=tilecol.*repmat(normals,1,3);
                end
            end
            meshcol=vcol;
            name=p.Results.coords;
            
            if length(vcol)==length(X)
                colormethod='interp';
            else
                colormethod='flat';
            end
            edgecolor='k';
            if ~p.Results.edges
                edgecolor='none';
            end
            flipped=obj.flipped;
            if ~p.Results.drawflips
                flipped=false(size(flipped));
            end
            patch('faces',obj.cut_mesh.T(flipped,:),'vertices',X,'facecolor','yellow','edgecolor','k','linewidth',2);%,'VertexNormals',normals);
            patch('faces',obj.cut_mesh.T(~flipped,:),'vertices',X,'facecolor',colormethod,'FaceVertexCData',vcol,'facealpha',1,'edgealpha',p.Results.edgealpha,'edgecolor',edgecolor,'linewidth',0.1);
            hold on
            for i=2:size(tileAs,3)
                Xtemp=X*tileAs(:,:,i)';
                
                patch('faces',obj.cut_mesh.T,'vertices',Xtemp,'facecolor',colormethod,'FaceVertexCData',tilecol,'facealpha',1,'edgealpha',1,'edgecolor','none','linewidth',0.1);%,'VertexNormals',normals);
                patch('faces',obj.cut_mesh.T,'vertices',Xtemp*1.001,'facecolor','none','edgealpha',0.1,'edgecolor','none','linewidth',0.1);
                if p.Results.tileboundary
                    obj.boundary.draw(Xtemp);
                end
            end
            
            
            if p.Results.boundary
                
                obj.boundary.draw(X);
                
            end
            
            S=[];
            axis off
            axis equal
            
        end
    end
    
end

