function contract = load_pacific_wind_data(cfg)
%LOAD_PACIFIC_WIND_DATA Load, orient, and reconstruct the 17-by-17 grid.

    windFile = locate_file(cfg, cfg.paths.windFileName, ...
        {'NCEP_U.dat', 'ncepu.dat'});
    coordinateFile = locate_file(cfg, cfg.paths.coordinateFileName, ...
        {'NCEPlatlon.dat', 'nceplatlon.dat'});
    U = readmatrix(windFile, 'FileType', 'text');
    coordinates = readmatrix(coordinateFile, 'FileType', 'text');

    if isequal(size(U), [cfg.data.expectedStations, ...
            cfg.data.expectedTimePoints])
        U = U.';
    end
    if isequal(size(coordinates), [2, cfg.data.expectedStations])
        coordinates = coordinates.';
    end
    if ~isequal(size(U), [cfg.data.expectedTimePoints, ...
            cfg.data.expectedStations])
        error('NCEP_U.dat has size %d-by-%d; expected 480-by-289.', ...
            size(U, 1), size(U, 2));
    end
    if ~isequal(size(coordinates), [cfg.data.expectedStations, 2])
        error(['NCEPlatlon.dat has size %d-by-%d; expected ', ...
            '289-by-2.'], size(coordinates, 1), size(coordinates, 2));
    end
    if any(~isfinite(coordinates), 'all')
        error('The coordinate file contains missing or nonfinite values.');
    end

    firstCouldBeLatitude = all(abs(coordinates(:, 1)) <= 90);
    secondCouldBeLatitude = all(abs(coordinates(:, 2)) <= 90);
    if firstCouldBeLatitude && ~secondCouldBeLatitude
        latitude = coordinates(:, 1);
        longitude = coordinates(:, 2);
        coordinateOrder = "latitude_longitude";
    elseif secondCouldBeLatitude && ~firstCouldBeLatitude
        latitude = coordinates(:, 2);
        longitude = coordinates(:, 1);
        coordinateOrder = "longitude_latitude";
    else
        error(['Could not determine latitude/longitude column order ', ...
            'unambiguously from coordinate ranges.']);
    end

    roundingScale = 1e8;
    latitudeRounded = round(latitude * roundingScale) / roundingScale;
    longitudeRounded = round(longitude * roundingScale) / roundingScale;
    latitudeLevels = unique(latitudeRounded, 'sorted');
    longitudeLevels = unique(longitudeRounded, 'sorted');
    if numel(latitudeLevels) ~= cfg.data.expectedGridSide || ...
            numel(longitudeLevels) ~= cfg.data.expectedGridSide
        error(['Coordinate file does not define a 17-by-17 regular grid: ', ...
            '%d latitude levels and %d longitude levels found.'], ...
            numel(latitudeLevels), numel(longitudeLevels));
    end
    [~, gridRow] = ismember(latitudeRounded, latitudeLevels);
    [~, gridColumn] = ismember(longitudeRounded, longitudeLevels);
    occupancy = accumarray([gridRow, gridColumn], 1, ...
        [cfg.data.expectedGridSide, cfg.data.expectedGridSide]);
    if any(occupancy ~= 1, 'all')
        error('The coordinate combinations do not form a complete grid.');
    end

    meanLatitude = mean(latitude);
    meanLongitude = mean(longitude);
    kilometersPerDegree = 111.195;
    xKm = (longitude - meanLongitude) * kilometersPerDegree * ...
        cosd(meanLatitude);
    yKm = (latitude - meanLatitude) * kilometersPerDegree;
    % Some MATLAB releases/configurations return compose output in a cell
    % representation.  Force a string column because later table matching
    % intentionally uses string equality.
    stationID = string(compose('G%02d_%02d', gridRow, gridColumn));
    originalColumn = (1:cfg.data.expectedStations)';
    siteTable = table(stationID, originalColumn, gridRow, gridColumn, ...
        latitude, longitude, xKm, yKm);

    contract.U = U;
    contract.coordinatesRaw = coordinates;
    contract.latitudeLevels = latitudeLevels;
    contract.longitudeLevels = longitudeLevels;
    contract.siteTable = siteTable;
    contract.xyKm = [xKm, yKm];
    contract.timeIndex = (1:cfg.data.expectedTimePoints)';
    contract.relativeHours = cfg.data.samplingHours * ...
        (0:(cfg.data.expectedTimePoints - 1))';
    contract.coordinateOrder = coordinateOrder;
    contract.windFile = windFile;
    contract.coordinateFile = coordinateFile;
    contract.responseName = cfg.data.responseName;
    contract.responseUnits = cfg.data.responseUnits;
end

function filePath = locate_file(cfg, preferredName, alternatives)
    candidates = strings(0, 1);
    roots = [string(cfg.paths.dataDir); ...
        string(fullfile(cfg.suiteRoot, 'Data', 'Pacific')); ...
        string(fullfile(pwd, 'Data')); ...
        string(fullfile(pwd, 'Data', 'Pacific')); string(pwd)];
    names = unique([string(preferredName); string(alternatives(:))], ...
        'stable');
    for r = 1:numel(roots)
        for n = 1:numel(names)
            candidates(end + 1, 1) = fullfile(roots(r), names(n)); %#ok<AGROW>
        end
    end
    filePath = '';
    for j = 1:numel(candidates)
        if isfile(candidates(j))
            filePath = char(candidates(j));
            return;
        end
    end
    error(['Could not find %s. Put both Pacific .dat files in ', ...
        '%s or supply them under Data/Pacific.'], preferredName, ...
        cfg.paths.dataDir);
end
