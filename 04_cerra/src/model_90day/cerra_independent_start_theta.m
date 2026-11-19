function theta=cerra_independent_start_theta(baseTheta,cfg,startID,seed)
%CERRA_INDEPENDENT_START_THETA Reproducible external multistart location.

    if startID==1
        theta=baseTheta;theta.nu=cfg.model.nu;return;
    end
    z=cerra_pack_constrained_parameters(baseTheta,cfg,[],[]);
    rng(seed+104729*startID,'twister');
    z=z+cfg.stage2.singleStartOptimizer.randomStartScale*randn(size(z));
    theta=cerra_unpack_constrained_parameters(z,cfg,[],[]);
end
