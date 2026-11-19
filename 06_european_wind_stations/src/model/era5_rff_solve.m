function solution = era5_rff_solve(factor, rightHandSide)
%ERA5_RFF_SOLVE Apply (Phi W Phi* + sigma_n^2 I)^(-1).

    noise = factor.noiseVariance;
    projected = factor.Z' * rightHandSide;
    correction = factor.L' \ (factor.L \ projected);
    solution = rightHandSide / noise - ...
        factor.Z * correction / (noise ^ 2);
end

