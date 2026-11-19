function summary=run_ghcnh_audit_scaling()
%RUN_GHCNH_AUDIT_SCALING Compare audit costs at fixed parameters and geometry.

    C=ghcnh_trace_config(); [cfg,~,~,A,E]=ghcnh_expanded_setup();
    active=readtable(fullfile(C.dataProcessedDir, ...
        'active_europe_station_inventory.csv'),'TextType','string');
    xyKm=project_europe_km(active.latitude,active.longitude);
    order=nested_space_filling_order(xyKm,1000);
    xyKm=xyKm(order,:); center=mean(xyKm,1);
    domain=(xyKm-center)/cfg.spatial.coordinateScaleKm;
    master=adaptive_basis_bank(domain,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    source=load(fullfile(C.resultsDir,'fold_01','algorithm1', ...
        'continuation_task.mat'),'task');
    theta=source.task.fit.theta;
    sizes=[200 429 671 800 1000]; rows=cell(numel(sizes),1);
    scalingDir=fullfile(E.resultsDir,'audit_scaling');
    if ~isfolder(scalingDir), mkdir(scalingDir); end
    logFile=fullfile(scalingDir,'progress.log');

    for r=1:numel(sizes)
        n=sizes(r); xy=domain(1:n,:);
        ghcnh_log(logFile,'audit_scaling',sprintf( ...
            'n=%d exact audit starting',n));
        timer=tic;
        [jExact,aExact,~]=exact_rank_audit(theta,xy,cfg,master, ...
            C,A,logFile,sprintf('scaling_n%d_exact',n));
        exactSeconds=toc(timer);
        timer=tic;
        [jTrace,aTrace,checks]=trace_rank_audit(theta,xy,cfg,master, ...
            C,A,logFile,sprintf('scaling_n%d_trace',n));
        traceSeconds=toc(timer);
        mExact=jExact.^2; mTrace=jTrace.^2;
        rows{r}=table(n,exactSeconds,traceSeconds, ...
            exactSeconds/traceSeconds,mean(jExact==jTrace), ...
            sum(jExact==C.referenceJ),sum(jTrace==C.referenceJ), ...
            sum(mExact.^3)/(numel(mExact)*C.referenceJ^6), ...
            sum(mTrace.^3)/(numel(mTrace)*C.referenceJ^6), ...
            'VariableNames',{'nSites','exactAuditSeconds', ...
            'traceAuditSeconds','auditTimeRatio','sameRankFraction', ...
            'exactReferenceFrequencies','traceReferenceFrequencies', ...
            'exactCubicRankLoad','traceCubicRankLoad'});
        writetable(aExact,fullfile(scalingDir,sprintf( ...
            'n_%04d_exact_audit.csv',n)));
        writetable(aTrace,fullfile(scalingDir,sprintf( ...
            'n_%04d_trace_audit.csv',n)));
        writetable(checks,fullfile(scalingDir,sprintf( ...
            'n_%04d_trace_spot_checks.csv',n)));
        profiles=table((0:numel(jExact)-1)',jExact,jTrace,mExact,mTrace, ...
            'VariableNames',{'frequencyIndex','exactJ','traceJ', ...
            'exactM','traceM'});
        writetable(profiles,fullfile(scalingDir,sprintf( ...
            'n_%04d_rank_profiles.csv',n)));
        writetable(vertcat(rows{1:r}),fullfile(scalingDir,'summary.csv'));
        ghcnh_log(logFile,'audit_scaling',sprintf( ...
            'n=%d done: exact %.2fs trace %.2fs same-rank %.3f', ...
            n,exactSeconds,traceSeconds,mean(jExact==jTrace)));
    end
    summary=vertcat(rows{:});
    save(fullfile(scalingDir,'audit_scaling.mat'), ...
        'summary','theta','sizes','order','-v7.3');
end

function order=nested_space_filling_order(xy,nTarget)
    z=rescale(xy,0,1); n=size(z,1); chosen=false(n,1);
    [~,first]=min(sum((z-0.5).^2,2)); chosen(first)=true;
    order=zeros(nTarget,1); order(1)=first;
    distance=sum((z-z(first,:)).^2,2);
    for j=2:nTarget
        distance(chosen)=-Inf; [~,next]=max(distance);
        chosen(next)=true; order(j)=next;
        distance=min(distance,sum((z-z(next,:)).^2,2));
    end
end

function xyKm=project_europe_km(latitude,longitude)
    if exist('projcrs','file')==2
        crs=projcrs(3035); [x,y]=projfwd(crs,latitude,longitude);
        xyKm=[x(:),y(:)]/1000; return
    end
    radius=6371.0088; lat=deg2rad(latitude(:)); lon=deg2rad(longitude(:));
    lat0=deg2rad(52); lon0=deg2rad(10); delta=lon-lon0;
    scale=sqrt(2./(1+sin(lat0).*sin(lat)+cos(lat0).*cos(lat).*cos(delta)));
    x=radius*scale.*cos(lat).*sin(delta);
    y=radius*scale.*(cos(lat0).*sin(lat)-sin(lat0).*cos(lat).*cos(delta));
    xyKm=[x,y];
end
