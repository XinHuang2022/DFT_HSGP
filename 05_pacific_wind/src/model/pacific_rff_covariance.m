function K = pacific_rff_covariance(PhiA, PhiB, weight, ...
        noiseVariance, sameLocations)
%PACIFIC_RFF_COVARIANCE Materialize an RFF covariance for diagnostics only.

    if nargin < 4 || isempty(noiseVariance)
        noiseVariance = 0;
    end
    if nargin < 5
        sameLocations = false;
    end
    K = PhiA * (weight(:) .* PhiB');
    K = real(K);
    if sameLocations
        K = (K + K') / 2 + noiseVariance * eye(size(K, 1));
    end
end
