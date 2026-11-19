function [L, jitterUsed] = eq9_chol_with_jitter(A)
%EQ9_CHOL_WITH_JITTER Cholesky factor with a small audited stabilization.

    A = real((A + A') / 2);
    scale = max(1, mean(abs(diag(A))));
    candidates = scale * [0, 1e-12, 1e-10, 1e-8, 1e-6];
    L = [];
    jitterUsed = NaN;
    for value = candidates
        [candidate, flag] = chol(A + value * eye(size(A)), 'lower');
        if flag == 0
            L = candidate;
            jitterUsed = value;
            return;
        end
    end
    error('Cholesky factorization failed after jitter escalation.');
end


