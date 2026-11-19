function K=cov_ref(theta,A,B,lambda,R,sameSites)
%COV_REF Numerical-Hankel half-spectrum plus optional white nugget.
    if nargin<6, sameSites=false; end
    g=(1-theta(2)^2)/(1+theta(2)^2-2*theta(2)*cos(lambda));
    kap=theta(3)*sqrt(1+theta(4)*sin(lambda/2)^2);
    K=theta(1)*g*rcorr(kap*distmat(A,B),R);
    if sameSites
        K=real((K+K')/2)+theta(5)*eye(size(A,1));
    end
end

