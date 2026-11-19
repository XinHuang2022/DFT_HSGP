function contract = load_cerra_wind_data(cfg)
%LOAD_CERRA_WIND_DATA Load native marine sites and the configured time window.

    windPath=fullfile(cfg.paths.dataDir,cfg.data.windFileName);
    sitePath=fullfile(cfg.paths.dataDir,cfg.data.siteFileName);
    timePath=fullfile(cfg.paths.dataDir,cfg.data.timeFileName);
    required={windPath,sitePath,timePath};
    for j=1:numel(required)
        if ~isfile(required{j}),error('Missing CERRA input: %s',required{j});end
    end
    windTable=readtable(windPath,'VariableNamingRule','preserve');
    site=readtable(sitePath,'VariableNamingRule','preserve');
    time=readtable(timePath,'VariableNamingRule','preserve');
    Uall=windTable{:,2:end};
    if size(Uall,1)~=cfg.data.expectedSourceTimePoints || ...
            size(Uall,2)~=cfg.data.expectedFullStations || ...
            any(~isfinite(Uall),'all')
        error('Full CERRA matrix must be finite and %d-by-%d.', ...
            cfg.data.expectedSourceTimePoints,cfg.data.expectedFullStations);
    end
    if height(site)~=cfg.data.expectedFullStations
        error('Full CERRA site manifest must contain %d sites.', ...
            cfg.data.expectedFullStations);
    end
    expectedHeaders=string(site.site_id(:));
    actualHeaders=string(windTable.Properties.VariableNames(2:end)).';
    if ~isequal(actualHeaders,expectedHeaders)
        error('Wind columns do not match site_manifest_full_native.csv.');
    end
    timeText=erase(string(time.time_utc),"Z");
    timeUtc=datetime(timeText,'InputFormat','yyyy-MM-dd''T''HH:mm:ss', ...
        'TimeZone','UTC');
    if numel(timeUtc)~=size(Uall,1)
        error('Time manifest has %d rows but the wind matrix has %d rows.', ...
            numel(timeUtc),size(Uall,1));
    end
    startUtc=datetime(erase(cfg.data.periodStartUtc,"Z"), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC');
    endUtc=datetime(erase(cfg.data.periodEndUtc,"Z"), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC');
    keepTime=timeUtc>=startUtc & timeUtc<=endUtc;
    if sum(keepTime)~=cfg.data.expectedTimePoints
        error('Configured period returned %d rather than %d time points.', ...
            sum(keepTime),cfg.data.expectedTimePoints);
    end
    if any(~isfinite(site.latitude)) || any(~isfinite(site.longitude)) || ...
            any(~isfinite(site.land_sea_mask)) || any(site.land_sea_mask>=0.5)
        error('Full-grid site manifest contains invalid or non-marine sites.');
    end
    lat0=mean(site.latitude);lon0=mean(site.longitude);
    earthRadiusKm=6371.0088;
    xKm=earthRadiusKm*cosd(lat0)*deg2rad(site.longitude-lon0);
    yKm=earthRadiusKm*deg2rad(site.latitude-lat0);
    siteTable=table(string(site.site_id),(1:height(site)).', ...
        site.native_y_index_zero_based+1,site.native_x_index_zero_based+1, ...
        site.latitude,site.longitude,xKm,yKm,site.land_sea_mask, ...
        'VariableNames',{'stationID','originalColumn','nativeRow','nativeColumn', ...
        'latitude','longitude','xKm','yKm','landSeaMask'});
    contract.U=Uall(keepTime,:);
    contract.siteTable=siteTable;
    contract.xyKm=[xKm,yKm];
    contract.timeUtc=timeUtc(keepTime);
    contract.timeIndex=(1:sum(keepTime)).';
    contract.relativeHours=hours(contract.timeUtc-contract.timeUtc(1));
    contract.sourceFiles=string(required(:));
    contract.responseName=string(cfg.data.responseName);
    contract.responseUnits=string(cfg.data.responseUnits);
end
