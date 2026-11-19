function D=distmat(A,B)
%DISTMAT Euclidean distance matrix without toolbox dependencies.
    D2=sum(A.^2,2)+sum(B.^2,2)'-2*(A*B');
    D=sqrt(max(0,D2));
end

