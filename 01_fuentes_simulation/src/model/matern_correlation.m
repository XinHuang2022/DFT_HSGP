function value = matern_correlation(distanceTimesKappa, eta)
%MATERN_CORRELATION Normalized Matérn correlation M_eta(x).

    x = max(0, distanceTimesKappa);
    if abs(eta - 0.5) < 1e-12
        value = exp(-x);
    elseif abs(eta - 1.5) < 1e-12
        value = (1 + x) .* exp(-x);
    elseif abs(eta - 2.5) < 1e-12
        value = (1 + x + x.^2 / 3) .* exp(-x);
    else
        value = ones(size(x));
        positive = x > 0;
        xp = x(positive);
        value(positive) = 2^(1 - eta) / gamma(eta) .* ...
            xp.^eta .* besselk(eta, xp);
    end
    value(x == 0) = 1;
end

