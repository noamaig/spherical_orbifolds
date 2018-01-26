function [ n ] = findHemisphere( p )
%find a hemisphere containing all points in p.
%output: the normal to the great circle of that hemisphere

%Las Vegas way - iterate over random normals until one is found

%a few canonical good guesses for the normal
N=[(mean(p)/norm(mean(p)));1 0 0;0 1 0;0 0 1;];

for iter=1:10000
    %first try the canonical guesses
    if ~isempty(N)
        n=N(1,:);
        N(1,:)=[];
    else
        %if trying just random,
        if mod(iter,2)==0
            n=normr(rand(1,3)-0.5);
        else
            inds=ceil(rand(2,1)*length(p));
            n=cross(p(inds(1),:),p(inds(2),:));
        end
    end
    if isOk(n)
        return;
    end
    n=-n;
    if isOk(n)
        return;
    end
    
    
end
error('points do not seem to lie all in same half-plane');
%old method - determinstic but some vertices lie on plane and need to fix
%that

    function ok=isOk(n)
        a=(p*n');
        tol=0;%1e-6;
        ok=all(a>=tol);
    end
end

