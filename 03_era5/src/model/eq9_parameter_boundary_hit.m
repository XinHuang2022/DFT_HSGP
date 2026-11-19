function [hit, minimumDistance] = eq9_parameter_boundary_hit(theta, cfg)
%EQ9_PARAMETER_BOUNDARY_HIT Natural-coordinate proximity to fitted bounds.

    names = {'signalVariance', 'alpha', 'beta', 'noiseVariance'};
    distances = zeros(numel(names) + 1, 1);
    for j = 1:numel(names)
        name = names{j};
        bounds = cfg.bounds.(name);
        fraction = (log(theta.(name)) - log(bounds(1))) / ...
            (log(bounds(2)) - log(bounds(1)));
        distances(j) = min(fraction, 1 - fraction);
    end
    delta = -0.5 * log(theta.epsilon);
    bounds = cfg.bounds.logSpatialScaleContrast;
    fraction = (delta - bounds(1)) / (bounds(2) - bounds(1));
    distances(end) = min(fraction, 1 - fraction);
    minimumDistance = min(distances);
    hit = minimumDistance < 0.01;
end

