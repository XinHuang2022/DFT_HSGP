function simulation = simulate_eq9_circulant_components( ...
    theta, xy, T, cfg, seed)
%SIMULATE_EQ9_CIRCULANT_COMPONENTS Circular latent process plus white noise.
%
% This SYN-4 simulator retains latent and noisy fields separately so that
% reconstruction and uncertainty can be scored for both scientific targets.

    rng(seed, 'twister');
    lambda = eq9_fourier_frequencies(T);
    [Fsignal, spectralDetails] = eq9_frequency_covariances( ...
        theta, xy, lambda, cfg, false);
    m = size(xy, 1);
    Dlatent = complex(zeros(T, m));
    jitters = zeros(T, 1);

    [L, jitters(1)] = eq9_chol_with_jitter(Fsignal(:, :, 1));
    Dlatent(1, :) = (L * randn(m, 1)).';
    if mod(T, 2) == 0
        lastInterior = T / 2;
        nyquistIndex = T / 2 + 1;
    else
        lastInterior = (T + 1) / 2;
        nyquistIndex = [];
    end
    for k = 2:lastInterior
        [L, jitters(k)] = eq9_chol_with_jitter(Fsignal(:, :, k));
        coefficient = L * (randn(m, 1) + 1i * randn(m, 1)) / sqrt(2);
        partner = T - k + 2;
        Dlatent(k, :) = coefficient.';
        Dlatent(partner, :) = conj(coefficient).';
        jitters(partner) = jitters(k);
    end
    if ~isempty(nyquistIndex)
        [L, jitters(nyquistIndex)] = eq9_chol_with_jitter( ...
            Fsignal(:, :, nyquistIndex));
        Dlatent(nyquistIndex, :) = (L * randn(m, 1)).';
    end

    complexLatent = ifft(Dlatent, [], 1) * sqrt(T);
    latent = real(complexLatent);
    noise = sqrt(theta.noiseVariance) * randn(T, m);
    Y = latent + noise;
    Dnoise = fft(noise, [], 1) / sqrt(T);

    simulation.Y = Y;
    simulation.latent = latent;
    simulation.noise = noise;
    simulation.D = Dlatent + Dnoise;
    simulation.Dlatent = Dlatent;
    simulation.Dnoise = Dnoise;
    simulation.Fsignal = Fsignal;
    simulation.lambda = lambda;
    simulation.theta = theta;
    simulation.xy = xy;
    simulation.seed = seed;
    simulation.choleskyJitter = jitters;
    simulation.maxImaginaryPart = max(abs(imag(complexLatent)), [], 'all');
    simulation.spectralDetails = spectralDetails;
end
