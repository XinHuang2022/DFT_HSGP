function solution = era5_hsgp_solve(factor, rightHandSide)
%ERA5_HSGP_SOLVE Apply a factor returned by era5_hsgp_factor.

    if factor.mode == "feature_qr"
        solution = era5_rff_solve(factor.core, rightHandSide);
    elseif factor.mode == "site_cholesky"
        solution = factor.L' \ (factor.L \ rightHandSide);
    else
        error('Unknown ERA5 HSGP factorization mode.');
    end
end

