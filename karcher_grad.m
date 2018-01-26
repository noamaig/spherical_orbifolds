function [gradu1,gradu2,gradu3,gradv1,gradv2,gradv3,o]=karcher_grad(u1,u2,u3,v1,v2,v3,w)

cross_a=u1.*v2 - u2.*v1;
cross_b=u1.*v3 - u3.*v1;
cross_c=u2.*v3 - u3.*v2;
dot_a=u1.*v1;
dot_b=u2.*v2;
dot_c=u3.*v3;
dot_sum=(dot_a + dot_b + dot_c);
cross_norm_squared=cross_a.^2 + cross_b.^2 + cross_c.^2;
cross_norm=sqrt(cross_norm_squared);
d=atan2(cross_norm,dot_sum);
dot_sum_squared=dot_sum.^2;
% gradu1 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((v1.*(cn).^(1./2))./(da + db + dc).^2 - (2.*v2.*(ca) + 2.*v3.*(cb))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);
%  
% gradu2 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((v2.*(cn).^(1./2))./(da + db + dc).^2 + (2.*v1.*(ca) - 2.*v3.*(cc))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);
%  
% gradu3 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((v3.*(cn).^(1./2))./(da + db + dc).^2 + (2.*v1.*(cb) + 2.*v2.*(cc))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);
%  
% gradv1 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((u1.*(cn).^(1./2))./(da + db + dc).^2 + (2.*u2.*(ca) + 2.*u3.*(cb))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);
%  
% gradv2 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((u2.*(cn).^(1./2))./(da + db + dc).^2 - (2.*u1.*(ca) - 2.*u3.*(cc))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);
%  
% gradv3 = -w.*(2.*atan((cn).^(1./2)./(da + db + dc)).*((u3.*(cn).^(1./2))./(da + db + dc).^2 - (2.*u1.*(cb) + 2.*u2.*(cc))./(2.*(da + db + dc).*(cn).^(1./2))))./((cn)./(da + db + dc).^2 + 1);

two_dot_sum_cross_norm=(2.*dot_sum.*cross_norm);
cross_norm_squared_over_dot_sum_squared_plus_one=(cross_norm_squared./dot_sum_squared + 1);
gradu1 = -w.*(2.*d.*((v1.*cross_norm)./dot_sum_squared - (2.*v2.*cross_a + 2.*v3.*cross_b)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;
 
gradu2 = -w.*(2.*d.*((v2.*cross_norm)./dot_sum_squared + (2.*v1.*cross_a - 2.*v3.*cross_c)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;
 
gradu3 = -w.*(2.*d.*((v3.*cross_norm)./dot_sum_squared + (2.*v1.*cross_b + 2.*v2.*cross_c)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;
 
gradv1 = -w.*(2.*d.*((u1.*cross_norm)./dot_sum_squared + (2.*u2.*cross_a + 2.*u3.*cross_b)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;
 
gradv2 = -w.*(2.*d.*((u2.*cross_norm)./dot_sum_squared - (2.*u1.*cross_a - 2.*u3.*cross_c)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;
 
gradv3 = -w.*(2.*d.*((u3.*cross_norm)./dot_sum_squared - (2.*u1.*cross_b + 2.*u2.*cross_c)./two_dot_sum_cross_norm))./cross_norm_squared_over_dot_sum_squared_plus_one;


o=sum(w.*d.^2);
%           gradu1 = -(2*atan(((u1*v2 - u2*v1)^2 + (u1*v3 - u3*v1)^2 + (u2*v3 - u3*v2)^2)^(1/2) /
%           (u1*v1 + u2*v2 + u3*v3))*(v1*u2^2 - u1*v2*u2 + v1*u3^2 - u1*v3*u3))/(((u1*v2 - u2*v1)^2 + (u1*v3 - u3*v1)^2 + (u2*v3 - u3*v2)^2)^(1/2)*(u1^2 + u2^2 + u3^2))


end
