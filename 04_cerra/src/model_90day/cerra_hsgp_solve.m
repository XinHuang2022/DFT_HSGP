function solution = cerra_hsgp_solve(factor, rightHandSide)
%CERRA_HSGP_SOLVE Apply a factor returned by cerra_hsgp_factor.

    if factor.mode == "feature_qr"
        solution = cerra_rff_solve(factor.core, rightHandSide);
    elseif factor.mode == "site_cholesky"
        solution = factor.L' \ (factor.L \ rightHandSide);
    else
        error('Unknown CERRA HSGP factorization mode.');
    end
end


