function result = run_ghcnh_metadata_audit()
%RUN_GHCNH_METADATA_AUDIT Screen European GHCNh stations without responses.

    C = ghcnh_trace_config();
    if ~isfolder(C.dataProcessedDir), mkdir(C.dataProcessedDir); end
    if ~isfolder(C.resultsDir), mkdir(C.resultsDir); end
    ghcnh_log(C.progressLog,'metadata','starting fixed-width station audit');
    stations = read_ghcnh_station_list(C.stationListFile);
    ghcnh_log(C.progressLog,'metadata',sprintf('parsed %d geolocated stations',height(stations)));

    inDomain = stations.latitude>=C.latitudeLimits(1) & ...
        stations.latitude<=C.latitudeLimits(2) & ...
        stations.longitude>=C.longitudeLimits(1) & ...
        stations.longitude<=C.longitudeLimits(2);
    europe = stations(inDomain,:);
    ghcnh_log(C.progressLog,'metadata',sprintf( ...
        '%d stations lie in lat [%g,%g], lon [%g,%g]',height(europe), ...
        C.latitudeLimits(1),C.latitudeLimits(2), ...
        C.longitudeLimits(1),C.longitudeLimits(2)));

    if isfile(C.inventoryFile)
        inventory = read_ghcnh_inventory_year(C.inventoryFile,C.year,C.progressLog);
        active = innerjoin(europe,inventory,'Keys','station');
        active.summerRecordCount = active.jun+active.jul+active.aug;
        active = active(active.summerRecordCount>=C.minimumSummerRecords,:);
        ghcnh_log(C.progressLog,'metadata',sprintf( ...
            '%d European stations have at least %d Jun-Aug %d records', ...
            height(active),C.minimumSummerRecords,C.year));
    else
        active=europe;
        active.summerRecordCount=nan(height(active),1);
        ghcnh_log(C.progressLog,'metadata', ...
            'monthly inventory absent; candidate selection uses all metadata stations');
    end

    % Deterministic space-filling ordering. The first pass selects at most
    % one station per 0.5-degree cell before filling from remaining stations.
    cellSize = 0.5;
    gx = floor((active.longitude-C.longitudeLimits(1))/cellSize);
    gy = floor((active.latitude-C.latitudeLimits(1))/cellSize);
    cellId = gy*1000+gx;
    [~,order] = sortrows([cellId, (1:height(active))']);
    ordered = active(order,:);
    orderedCell = cellId(order);
    [~,first] = unique(orderedCell,'stable');
    remaining = setdiff((1:height(ordered))',first,'stable');
    candidateOrder = [first; remaining];
    nCandidate = min(C.downloadCandidateCount,height(ordered));
    candidates = ordered(candidateOrder(1:nCandidate),:);
    candidates.downloadOrder = (1:height(candidates))';

    writetable(europe,fullfile(C.dataProcessedDir,'europe_station_inventory.csv'));
    writetable(active,fullfile(C.dataProcessedDir,'active_europe_station_inventory.csv'));
    writetable(candidates,fullfile(C.dataProcessedDir,'download_candidates.csv'));
    save(fullfile(C.dataProcessedDir,'metadata_audit.mat'),'C','stations', ...
        'europe','active','candidates');

    byCountry = groupsummary(europe,'country');
    byCountry = sortrows(byCountry,'GroupCount','descend');
    writetable(byCountry,fullfile(C.dataProcessedDir,'europe_counts_by_country.csv'));
    ghcnh_log(C.progressLog,'metadata',sprintf( ...
        'saved %d deterministic download candidates',height(candidates)));
    result = struct('allStations',height(stations),'europeStations', ...
        height(europe),'activeEuropeStations',height(active), ...
        'downloadCandidates',height(candidates), ...
        'candidateFile',fullfile(C.dataProcessedDir,'download_candidates.csv'));
end
