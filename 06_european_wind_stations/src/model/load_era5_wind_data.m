function contract = load_era5_wind_data(cfg)
%LOAD_ERA5_WIND_DATA Read the fixed 17-by-17 North Sea ERA5 data contract.

    windPath = fullfile(cfg.paths.dataDir, cfg.data.windFileName);
    sitePath = fullfile(cfg.paths.dataDir, cfg.data.siteFileName);
    timePath = fullfile(cfg.paths.dataDir, cfg.data.timeFileName);
    required = {windPath, sitePath, timePath};
    for j = 1:numel(required)
        if ~isfile(required{j}), error('Missing ERA5 input: %s', required{j}); end
    end
    windTable = readtable(windPath, 'VariableNamingRule', 'preserve');
    site = readtable(sitePath, 'VariableNamingRule', 'preserve');
    time = readtable(timePath, 'VariableNamingRule', 'preserve');
    U = windTable{:, 2:end};
    if ~isequal(size(U), [cfg.data.expectedTimePoints, ...
            cfg.data.expectedStations]) || any(~isfinite(U), 'all')
        error('ERA5 u10 matrix must be finite and 360-by-289.');
    end
    expectedHeaders = string(site.site_id(:));
    actualHeaders = string(windTable.Properties.VariableNames(2:end)).';
    if ~isequal(actualHeaders, expectedHeaders)
        error('Wind columns do not match site_manifest_17x17.csv.');
    end
    if height(site) ~= cfg.data.expectedStations || ...
            numel(unique(site.row)) ~= cfg.data.expectedGridSide || ...
            numel(unique(site.column)) ~= cfg.data.expectedGridSide
        error('ERA5 site manifest is not the expected 17-by-17 grid.');
    end
    if any(site.land_sea_mask ~= 0)
        error('This experiment requires the fixed all-ocean North Sea grid.');
    end
    lat0 = mean(site.latitude); lon0 = mean(site.longitude);
    earthRadiusKm = 6371.0088;
    xKm = earthRadiusKm * cosd(lat0) * deg2rad(site.longitude - lon0);
    yKm = earthRadiusKm * deg2rad(site.latitude - lat0);
    siteTable = table(string(site.site_id), (1:height(site)).', ...
        site.row, site.column, site.latitude, site.longitude, ...
        xKm, yKm, site.land_sea_mask, 'VariableNames', ...
        {'stationID','originalColumn','gridRow','gridColumn', ...
         'latitude','longitude','xKm','yKm','landSeaMask'});
    timeText = erase(string(time.time_utc), "Z");
    timeUtc = datetime(timeText, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'UTC');
    if any(isnat(timeUtc)) || height(time) ~= size(U, 1)
        error('Invalid ERA5 time manifest.');
    end
    contract.U = U;
    contract.siteTable = siteTable;
    contract.xyKm = [xKm, yKm];
    contract.timeUtc = timeUtc;
    contract.timeIndex = (1:size(U, 1)).';
    contract.relativeHours = hours(timeUtc - timeUtc(1));
    contract.sourceFiles = string(required(:));
    contract.responseName = string(cfg.data.responseName);
    contract.responseUnits = string(cfg.data.responseUnits);
end
