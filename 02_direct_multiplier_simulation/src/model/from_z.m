function theta=from_z(z,C)
%FROM_Z Logistic map from optimizer coordinates to parameter bounds.
    p=1./(1+exp(-z(:)'));
    theta=C.lower+(C.upper-C.lower).*p;
end

