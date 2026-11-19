function subset=prepare_ghcnh_expanded_panel()
%PREPARE_GHCNH_EXPANDED_PANEL Build the strict complete-case expanded panel.

    C=ghcnh_trace_config();
    files=dir(fullfile(C.dataProcessedDir,'wind_panel_expanded_*.mat'));
    if isempty(files), error('Expanded completeness output is absent.'); end
    [~,latest]=max([files.datenum]);
    S=load(fullfile(files(latest).folder,files(latest).name), ...
        'coverage','Y','targetTimes');
    use=find(S.coverage.complete);
    siteTable=S.coverage(use,:); U=S.Y(use,:).';
    xyKm=project_europe_km(siteTable.latitude,siteTable.longitude);

    subset=struct(); subset.U=U; subset.xyKm=xyKm;
    subset.siteTable=siteTable;
    subset.siteTable.stationID=string(siteTable.station);
    subset.siteTable.gridRow=(1:height(siteTable))';
    subset.siteTable.gridColumn=zeros(height(siteTable),1);
    subset.targetTimes=S.targetTimes;
    subset.timeIndex=(0:size(U,1)-1)';
    subset.relativeHours=hours(S.targetTimes-S.targetTimes(1));
    subset.responseName='wind_speed'; subset.responseUnits='m/s';
    subset.dataSource='NOAA GHCNh expanded strict complete-case panel';
    output=fullfile(C.dataProcessedDir,'ghcnh_complete_panel_expanded.mat');
    save(output,'C','subset','-v7.3');
    writetable(subset.siteTable,fullfile(C.dataProcessedDir, ...
        'ghcnh_complete_station_manifest_expanded.csv'));
    ghcnh_log(C.progressLog,'expanded_prepare',sprintf( ...
        'prepared T=%d by n=%d strict complete-case panel', ...
        size(U,1),size(U,2)));
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
