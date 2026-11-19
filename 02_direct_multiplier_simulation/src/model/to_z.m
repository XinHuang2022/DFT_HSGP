function z=to_z(theta,C)
%TO_Z Map bounded scientific parameters to an unconstrained vector.
    p=(theta-C.lower)./(C.upper-C.lower);
    p=min(1-1e-10,max(1e-10,p));
    z=log(p./(1-p)); z=z(:);
end

