    
        function [O,G]=objective_and_grad(obj,X)
            assert(size(X,2)==1);
            assert(~any(isnan(X(:))),'initial guess has NaN''s!');
            %given a vector s, return the jacboain G of the karcher
            %objectives, that is
            %G_ij is the hyperbolic distance between the vertices of the
            %i'th edge, derived according to the j'th variable
            %karcher
            
            
            
            I=obj.I;
            J=obj.J;
            
            u1=X(I*3-2);
            u2=X(I*3-1);
            u3=X(I*3);
            v1=X(J*3-2);
            v2=X(J*3-1);
            v3=X(J*3);
            [gradu1,gradu2,gradu3,gradv1,gradv2,gradv3,O]=karcher_grad(u1,u2,u3,v1,v2,v3,obj.W);
            %%
            ff = find(isnan(gradu1));
            if ~isempty(ff)
                %warning('nans in gradient!');
                gradu1(ff) = 0;
                gradu2(ff) = 0;
                gradu3(ff) = 0;
                gradv1(ff) = 0;
                gradv2(ff) = 0;
                gradv3(ff) = 0;
            end
%            	%% sanity check!!
%             gg=[gradu1 gradu2 gradu3];
%             xx=[u1 u2 u3];
%             aa=dot(normr(gg)',normr(xx)');
%             aa=max(abs(aa));
%             if(aa)>1e-6
%                 aa
%             end
            %%
            w=obj.W;
            
            Ju1=I*3-2;
            Ju2=I*3-1;
            Ju3=I*3;
            Jv1=J*3-2;
            Jv2=J*3-1;
            Jv3=J*3;
            
            J_k=[Ju1;Ju2;Ju3;Jv1;Jv2;Jv3];
            
            V_k=[gradu1;gradu2;gradu3;gradv1;gradv2;gradv3];
            
            
            if max(J_k)<length(obj.adj)
                J_k=[J_k;length(obj.adj)];
                V_k=[V_k;0];
            end
            
            G = accumarray(J_k, V_k);
            
%             inds=find(obj.cones);
%             G(inds*2-1)=0;
%             G(inds*2)=0;
            %make sure the "redundant" vertices (cones + right side) have
            %no effect on jacobian.
            %             rs=find(obj.right_side_vertices);
            %             cones=find(obj.cones);
            %             assert(all(all(G(rs*2-1)==0)));
            %             assert(all(all(G(rs*2)==0)));
            %             assert(all(all(G(cones*2-1)==0)));
            %             assert(all(all(G(cones*2)==0)));
        end
    