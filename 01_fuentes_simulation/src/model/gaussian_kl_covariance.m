function [klValue, jitterP, jitterQ] = gaussian_kl_covariance(P, Q, factor)
%GAUSSIAN_KL_COVARIANCE KL{N(0,P) || N(0,Q)} for real/complex blocks.
%
% factor=1/2 gives the real-Gaussian expression (zero and Nyquist DFT
% frequencies). factor=1 gives the proper complex-Gaussian expression used
% at interior positive frequencies.

    if nargin < 3
        factor = 0.5;
    end
    P = real((P + P') / 2);
    Q = real((Q + Q') / 2);
    [LP, jitterP] = stable_chol(P);
    [LQ, jitterQ] = stable_chol(Q);
    n = size(P, 1);
    logdetP = 2 * sum(log(diag(LP)));
    logdetQ = 2 * sum(log(diag(LQ)));
    traceTerm = trace(LQ' \ (LQ \ P));
    klValue = factor * real(traceTerm - n + logdetQ - logdetP);
    if klValue < 0 && klValue > -1e-9
        klValue = 0;
    elseif klValue < -1e-9
        error('Computed a materially negative covariance KL divergence.');
    end
end

function [L, jitter] = stable_chol(A)
    scale = max(1, mean(diag(A)));
    jitter = 0;
    for attempt = 1:10
        [L, flag] = chol(A + jitter * eye(size(A)), 'lower');
        if flag == 0
            return;
        end
        if jitter == 0
            jitter = 1e-12 * scale;
        else
            jitter = 10 * jitter;
        end
    end
    error('Cholesky factorization failed after adaptive jitter.');
end

