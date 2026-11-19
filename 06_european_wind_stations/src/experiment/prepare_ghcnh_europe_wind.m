function subset = prepare_ghcnh_europe_wind(maxStations)
%PREPARE_GHCNH_EUROPE_WIND Create the complete-record European wind panel.

    C=ghcnh_trace_config();
    if nargin<1 || isempty(maxStations), maxStations=C.downloadCandidateCount; end
    file=fullfile(C.dataProcessedDir,sprintf('wind_panel_%04d.mat',maxStations));
    if ~isfile(file), error('Completeness output is absent: %s',file); end
    S=load(file,'coverage','Y','targetTimes');
    use=find(S.coverage.complete);
    if numel(use)<C.targetStationRange(1)
        error('Only %d complete stations; at least %d were prescribed.', ...
            numel(use),C.targetStationRange(1));
    end
    if numel(use)>C.targetStationRange(2)
        use=space_filling_subset(S.coverage(use,:),C.targetStationRange(2));
    end
    siteTable=S.coverage(use,:);
    U=S.Y(use,:).';
    xyKm=project_europe_km(siteTable.latitude,siteTable.longitude);

    subset=struct();
    subset.U=U;
    subset.xyKm=xyKm;
    subset.siteTable=siteTable;
    subset.siteTable.stationID=string(siteTable.station);
    subset.siteTable.gridRow=(1:height(siteTable))';
    subset.siteTable.gridColumn=zeros(height(siteTable),1);
    subset.targetTimes=S.targetTimes;
    subset.timeIndex=(0:size(U,1)-1)';
    subset.relativeHours=hours(S.targetTimes-S.targetTimes(1));
    subset.responseName='wind_speed';
    subset.responseUnits='m/s';
    subset.dataSource='NOAA GHCNh';

    save(fullfile(C.dataProcessedDir,'ghcnh_complete_panel.mat'), ...
        'C','subset','-v7.3');
    writetable(subset.siteTable,fullfile(C.dataProcessedDir, ...
        'ghcnh_complete_station_manifest.csv'));
    ghcnh_log(C.progressLog,'prepare',sprintf( ...
        'prepared T=%d by n=%d complete wind panel',size(U,1),size(U,2)));
end

function selected=space_filling_subset(tbl,nTarget)
    x=rescale(tbl.longitude,0,1); y=rescale(tbl.latitude,0,1);
    chosen=false(height(tbl),1);
    [~,first]=min((x-0.5).^2+(y-0.5).^2); chosen(first)=true;
    distance=(x-x(first)).^2+(y-y(first)).^2;
    for j=2:nTarget
        distance(chosen)=-Inf;
        [~,next]=max(distance); chosen(next)=true;
        distance=max(0,min(distance,(x-x(next)).^2+(y-y(next)).^2));
    end
    selected=find(chosen);
end

function xyKm=project_europe_km(latitude,longitude)
    if exist('projcrs','file')==2
        crs=projcrs(3035);
        [x,y]=projfwd(crs,latitude,longitude);
        xyKm=[x(:),y(:)]/1000;
        return;
    end
    % Spherical Lambert azimuthal equal-area fallback, centered on Europe.
    radius=6371.0088; lat=deg2rad(latitude(:)); lon=deg2rad(longitude(:));
    lat0=deg2rad(52); lon0=deg2rad(10); delta=lon-lon0;
    scale=sqrt(2./(1+sin(lat0).*sin(lat)+cos(lat0).*cos(lat).*cos(delta)));
    x=radius*scale.*cos(lat).*sin(delta);
    y=radius*scale.*(cos(lat0).*sin(lat)-sin(lat0).*cos(lat).*cos(delta));
    xyKm=[x,y];
end
