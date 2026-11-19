function P=features(xy,B)
%FEATURES Evaluate normalized tensor-product Dirichlet eigenfunctions.
    q=xy-B.center;
    if any(abs(q)>B.L+1e-12,'all'), error('Site outside HSGP domain.'); end
    P=ones(size(xy,1),B.M);
    for d=1:2
        P=P.*sin((q(:,d)+B.L(d)).*B.omega(:,d)')./sqrt(B.L(d));
    end
end

