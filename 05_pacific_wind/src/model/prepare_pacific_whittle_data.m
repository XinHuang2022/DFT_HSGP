function prepared = prepare_pacific_whittle_data(Y, xyCovariance, design)
%PREPARE_PACIFIC_WHITTLE_DATA Cache data/design DFTs outside optimization.

    [T, nStations] = size(Y);
    if any(~isfinite(Y), 'all') || size(xyCovariance, 1) ~= nStations
        error('Pacific Whittle fitting requires complete T-by-site data.');
    end
    if design.T ~= T || design.nStations ~= nStations
        error('Mean design does not match the Pacific data dimensions.');
    end
    if mod(T, 2) == 0
        positiveIndices = 2:T/2;
        endpointIndices = [1, T/2 + 1];
    else
        positiveIndices = 2:(T + 1)/2;
        endpointIndices = 1;
    end
    prepared.D = fft(Y, [], 1) / sqrt(T);
    prepared.Afrequency = design.Afrequency;
    prepared.xy = xyCovariance;
    prepared.T = T;
    prepared.nStations = nStations;
    prepared.nCoefficients = design.nCoefficients;
    prepared.endpointIndices = endpointIndices;
    prepared.positiveIndices = positiveIndices;
    prepared.usedIndices = [endpointIndices, positiveIndices];
end
