function bank = make_pacific_rff_bank(thetaReference, nFeatures, cfg, seed)
%MAKE_PACIFIC_RFF_BANK Fixed paired spatial-frequency importance sample.

% q(w) is an equal-or-user-weighted mixture of normalized Matérn spectral
% densities at the low, transition, and high temporal-frequency spatial
% scales implied by thetaReference.  Frequencies occur as +/- pairs, making
% every approximated covariance real and positive semidefinite.

    if mod(nFeatures, 2) ~= 0 || nFeatures < 4
        error('The paired RFF feature count must be even and at least four.');
    end
    d = cfg.model.dimension;
    eta = thetaReference.nu - d / 2;
    if eta <= 0
        error('RFF Matérn eta must be positive.');
    end
    kappaLow = thetaReference.alpha;
    kappaMiddle = thetaReference.alpha * sqrt(2 / ...
        (1 + thetaReference.epsilon));
    kappaHigh = thetaReference.alpha / sqrt(thetaReference.epsilon);
    proposalKappa = [kappaLow, kappaMiddle, kappaHigh];
    mixture = cfg.rff.proposalMixtureWeights(:).';
    if numel(mixture) ~= 3 || any(mixture <= 0)
        error('RFF proposal mixture must contain three positive weights.');
    end
    mixture = mixture / sum(mixture);

    rng(seed, 'twister');
    nBase = nFeatures / 2;
    uniform = rand(nBase, 1);
    component = 1 + (uniform > mixture(1)) + ...
        (uniform > sum(mixture(1:2)));
    % p_kappa is a multivariate t density with 2*eta degrees of freedom
    % and scale kappa^2/(2*eta).  If V~chi2(2*eta), omega=kappa Z/sqrt(V).
    chiSquared = 2 * randg(eta, nBase, 1);
    baseOmega = randn(nBase, d) .* proposalKappa(component).' ./ ...
        sqrt(chiSquared);
    omega = [baseOmega; -baseOmega];
    component = [component; component];
    componentDensity = pacific_matern_spectral_density( ...
        omega, proposalKappa, eta);
    proposalDensity = componentDensity * mixture.';
    if any(~isfinite(proposalDensity)) || any(proposalDensity <= 0)
        error('RFF proposal density is not positive and finite.');
    end

    bank.version = "paired_matern_mixture_v1";
    bank.seed = seed;
    bank.nFeatures = nFeatures;
    bank.dimension = d;
    bank.eta = eta;
    bank.omega = omega;
    bank.proposalDensity = proposalDensity;
    bank.sampledComponent = component;
    bank.proposalKappa = proposalKappa;
    bank.mixtureWeights = mixture;
    bank.thetaReference = thetaReference;
end
