function candidates=prepare_ghcnh_expanded_manifest()
%PREPARE_GHCNH_EXPANDED_MANIFEST Preserve baseline order and add all candidates.

    C=ghcnh_trace_config();
    activeFile=fullfile(C.dataProcessedDir,'active_europe_station_inventory.csv');
    baselineFile=fullfile(C.dataProcessedDir,'download_candidates.csv');
    active=readtable(activeFile,'TextType','string');
    baseline=readtable(baselineFile,'TextType','string');

    cellSize=0.5;
    gx=floor((active.longitude-C.longitudeLimits(1))/cellSize);
    gy=floor((active.latitude-C.latitudeLimits(1))/cellSize);
    cellId=gy*1000+gx;
    [~,order]=sortrows([cellId,(1:height(active))']);
    ordered=active(order,:); orderedCell=cellId(order);
    [~,first]=unique(orderedCell,'stable');
    remaining=setdiff((1:height(ordered))',first,'stable');
    candidateOrder=[first;remaining];
    candidates=ordered(candidateOrder,:);
    candidates.downloadOrder=(1:height(candidates))';

    if height(baseline)>height(candidates) || ...
            any(baseline.station~=candidates.station(1:height(baseline)))
        error('ExpandedManifest:BaselineMismatch', ...
            'Expanded ordering does not preserve the original station manifest.');
    end
    outputFile=fullfile(C.dataProcessedDir,sprintf( ...
        'download_candidates_expanded_%04d.csv',height(candidates)));
    writetable(candidates,outputFile);
    ghcnh_log(C.progressLog,'expanded_manifest',sprintf( ...
        'saved %d candidates; original %d-station prefix verified', ...
        height(candidates),height(baseline)));
end
