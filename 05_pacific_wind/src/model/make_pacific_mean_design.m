function design = make_pacific_mean_design(T, xyStandardized)
%MAKE_PACIFIC_MEAN_DESIGN Spatial plane, constant through time.

    nStations = size(xyStandardized, 1);
    spatial = [ones(nStations, 1), xyStandardized];
    p = size(spatial, 2);
    B = zeros(T, nStations, p);
    for j = 1:p
        B(:, :, j) = repmat(spatial(:, j).', T, 1);
    end
    design.model = "P0_spatial_plane";
    design.B = B;
    design.Afrequency = fft(B, [], 1) / sqrt(T);
    design.coefficientNames = ["intercept"; "east"; "north"];
    design.nCoefficients = p;
    design.nStations = nStations;
    design.T = T;
end
