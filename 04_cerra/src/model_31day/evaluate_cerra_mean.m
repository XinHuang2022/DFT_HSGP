function fittedMean = evaluate_cerra_mean(design, eta)
%EVALUATE_CERRA_MEAN Evaluate the time-domain spatial-plane mean.

    eta = eta(:);
    if numel(eta) ~= design.nCoefficients
        error('Mean coefficient vector has incorrect length.');
    end
    fittedMean = zeros(design.T, design.nStations);
    for j = 1:design.nCoefficients
        fittedMean = fittedMean + design.B(:, :, j) * eta(j);
    end
end


