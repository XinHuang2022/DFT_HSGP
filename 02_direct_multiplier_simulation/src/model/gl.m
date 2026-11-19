function [x,w]=gl(n,a,b)
%GL Gauss-Legendre nodes and weights on [a,b] (Golub-Welsch).
    j=(1:n-1)'; beta=j./sqrt(4*j.^2-1);
    [V,L]=eig(diag(beta,1)+diag(beta,-1));
    [z,ix]=sort(diag(L)); V=V(:,ix);
    wg=2*(V(1,:)').^2;
    x=(a+b)/2+(b-a)*z/2;
    w=(b-a)*wg/2;
end

