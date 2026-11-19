function S=mult(theta,omega,lambda)
%MULT Spatial covariance-operator multiplier S_theta(omega,lambda).
% theta=[vf,rho,kappa0,gamma,vn]; omega can be a vector of radial norms.
    vf=theta(1); rho=theta(2); k0=theta(3); gam=theta(4);
    g=(1-rho^2)./(1+rho^2-2*rho*cos(lambda));
    kap=k0.*sqrt(1+gam*sin(lambda/2).^2);
    omega=omega(:); lambda=lambda(:).';
    kap=kap(:).'; g=g(:).';
    S=vf.*g.*(9*sqrt(3)./(kap.^2)).* ...
        (1+(omega./kap).^3).^(-2);
end

