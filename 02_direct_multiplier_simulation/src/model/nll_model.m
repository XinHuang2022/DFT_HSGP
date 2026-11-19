function [value,detail]=nll_model(z,Y,xy,C,R,engine,B,cache)
%NLL_MODEL One-sided exact circular Gaussian/Whittle objective.
    if nargin<8, cache=[]; end
    theta=from_z(z,C); [T,n]=size(Y); D=fft(Y,[],1)/sqrt(T);
    lam=freqs(T); value=0; maxJ=0;
    if mod(T,2)==0, ends=[1,T/2+1]; interior=2:T/2;
    else, ends=1; interior=2:(T+1)/2; end
    if strcmp(engine,'hsgp')
        if isempty(cache)
            P=features(xy,B); G=P'*P;
        else
            P=cache.P; G=cache.G;
        end
    end
    try
        for k=ends
            [F,L]=factor_at(k); d=real(D(k,:)).';
            value=value+.5*(F+real(d'*solve_at(L,d))+n*log(2*pi));
            maxJ=max(maxJ,L.jit);
        end
        for k=interior
            [F,L]=factor_at(k); d=D(k,:).';
            value=value+F+real(d'*solve_at(L,d))+n*log(pi);
            maxJ=max(maxJ,L.jit);
        end
        if ~isfinite(value), error('Nonfinite objective.'); end
        detail.theta=theta; detail.maxJitter=maxJ;
    catch ME
        value=realmax/1e100; detail.error=ME.message; detail.maxJitter=Inf;
    end

    function [ld,F]=factor_at(k)
        if strcmp(engine,'exact')
            K=cov_ref(theta,xy,xy,lam(k),R,true);
            [Lc,j]=cholj(K); F.kind='dense'; F.L=Lc; F.jit=j;
            ld=2*sum(log(diag(Lc)));
        else
            w=mult(theta,B.radial,lam(k));
            F=lrfactor(P,w,theta(5),G); ld=F.logdet;
        end
    end
    function x=solve_at(F,b)
        if strcmp(F.kind,'dense'), x=F.L'\(F.L\b); else, x=lrsolve(F,b); end
    end
end
