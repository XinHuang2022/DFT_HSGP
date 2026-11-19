function basis = make_cerra_hsgp_basis(domainXy,boundaryFactor,nFeatures,cfg)
%MAKE_CERRA_HSGP_BASIS Lowest radial Laplacian modes on a rectangle.

    if cfg.model.dimension~=2 || size(domainXy,2)~=2 || ...
            any(~isfinite(domainXy),'all')
        error('CERRA HSGP requires finite two-dimensional coordinates.');
    end
    if ~(isscalar(boundaryFactor)&&isfinite(boundaryFactor)&&boundaryFactor>1)
        error('boundaryFactor must be finite and larger than one.');
    end
    if ~(isscalar(nFeatures)&&nFeatures>=1&&fix(nFeatures)==nFeatures)
        error('nFeatures must be a positive integer.');
    end
    lower=min(domainXy,[],1);upper=max(domainXy,[],1);
    center=(lower+upper)/2;supportHalfWidth=(upper-lower)/2;
    if any(supportHalfWidth<=0),error('HSGP domain must span both coordinates.');end
    L=boundaryFactor*supportHalfWidth;
    geometricMean=sqrt(prod(L));
    maxIndex=max(2,ceil(2.5*sqrt(nFeatures)*L/geometricMean));
    [j1,j2]=ndgrid(1:maxIndex(1),1:maxIndex(2));
    candidate=[j1(:),j2(:)];
    candidateOmega=pi*candidate./(2*L);
    candidateLambda=sum(candidateOmega.^2,2);
    [~,order]=sortrows([candidateLambda,candidate(:,1),candidate(:,2)],[1 2 3]);
    if numel(order)<nFeatures,error('Internal radial-mode grid is too small.');end
    indices=candidate(order(1:nFeatures),:);
    omega=pi*indices./(2*L);
    eigenvalue=sum(omega.^2,2);
    basis.method='HSGP';
    basis.boundaryCondition='dirichlet';
    basis.truncation='radial_laplacian_eigenvalue';
    basis.dimension=2;
    basis.boundaryFactor=boundaryFactor;
    basis.nFeatures=nFeatures;
    basis.center=center;
    basis.supportHalfWidth=supportHalfWidth;
    basis.domainHalfWidth=L;
    basis.indices=indices;
    basis.maximumIndexByDimension=max(indices,[],1);
    basis.omegaByDimension=omega;
    basis.laplaceEigenvalue=eigenvalue;
    basis.radialFrequency=sqrt(eigenvalue);
    basis.maximumRadialFrequency=max(basis.radialFrequency);
end
