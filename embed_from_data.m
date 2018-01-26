function e= embed_from_data( Vo,To,cone_order,cone_vertex_indices)
%embed the given mesh (Vo,To) according to the cone order and the indices
%of the cone vertices
os=createOrbifoldStructure(Vo,To,cone_order,cone_vertex_indices);
[boundary,cutter]=os.generate();
e=embed_from_boundary(boundary,cutter);

end

