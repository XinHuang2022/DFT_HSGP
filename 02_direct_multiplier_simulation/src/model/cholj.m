function [L,jit]=cholj(K)
%CHOLJ Cholesky with recorded, scale-aware numerical jitter.
    K=real((K+K')/2); s=max(1,mean(abs(diag(K))));
    cand=s*[0,1e-12,1e-10,1e-8,1e-6];
    for jit=cand
        [L,p]=chol(K+jit*eye(size(K)),'lower');
        if p==0, return; end
    end
    error('Cholesky failed after jitter ladder.');
end
