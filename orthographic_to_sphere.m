function  V  = orthographic_to_sphere( v )
Z=real(sqrt(1-sum(v.^2,2)));

V=[v Z];

end

