function F=lrfactor(P,w,vn,G)
%LRFACTOR Cached-Gram Woodbury factor; dense fallback when M >= n.
    n=size(P,1); m=size(P,2); w=w(:);
    if m<n
        A=diag(1./w)+G/vn;
        [L,jit]=cholj(A);
        F.kind='woodbury'; F.P=P; F.L=L; F.vn=vn; F.jit=jit;
        F.logdet=n*log(vn)+sum(log(w))+2*sum(log(diag(L)));
    else
        K=(P.*w')*P'+vn*eye(n);
        [L,jit]=cholj(K);
        F.kind='dense'; F.L=L; F.jit=jit;
        F.logdet=2*sum(log(diag(L)));
    end
end

