function summary = download_ghcnh_candidates(maxStations)
%DOWNLOAD_GHCNH_CANDIDATES Download compact GHCNh station-year Parquet files.

    C = ghcnh_trace_config();
    if nargin<1 || isempty(maxStations), maxStations=C.downloadCandidateCount; end
    candidates = readtable(fullfile(C.dataProcessedDir,'download_candidates.csv'), ...
        'TextType','string');
    maxStations=min(maxStations,height(candidates));
    candidates=candidates(1:maxStations,:);
    if ~isfolder(C.stationDataDir), mkdir(C.stationDataDir); end
    if ~isfolder(C.resultsDir), mkdir(C.resultsDir); end
    ghcnh_log(C.progressLog,'download',sprintf( ...
        'starting/resuming %d station files with %d workers', ...
        maxStations,C.downloadWorkers));

    pool=gcp('nocreate');
    if isempty(pool) || pool.NumWorkers~=C.downloadWorkers
        if ~isempty(pool), delete(pool); end
        pool=parpool('local',C.downloadWorkers);
    end
    futures=parallel.FevalFuture.empty(maxStations,0);
    for i=1:maxStations
        station=char(candidates.station(i));
        url=sprintf(C.stationParquetUrl,C.year,station,C.year);
        destination=fullfile(C.stationDataDir, ...
            sprintf('GHCNh_%s_%d.parquet',station,C.year));
        futures(i)=parfeval(pool,@download_one,1,url,destination,station);
    end

    station=strings(maxStations,1); status=strings(maxStations,1);
    bytes=zeros(maxStations,1); elapsed=zeros(maxStations,1);
    for completed=1:maxStations
        [index,out]=fetchNext(futures);
        station(index)=string(out.station); status(index)=string(out.status);
        bytes(index)=out.bytes; elapsed(index)=out.elapsedSeconds;
        if mod(completed,10)==0 || completed==maxStations
            ghcnh_log(C.progressLog,'download',sprintf( ...
                '%d/%d complete; ok=%d, missing=%d, failed=%d; %.1f MB', ...
                completed,maxStations,sum(status=="downloaded" | status=="existing"), ...
                sum(status=="missing"),sum(status=="failed"),sum(bytes)/2^20));
        end
    end
    summary=table(station,status,bytes,elapsed, ...
        'VariableNames',{'station','status','bytes','elapsedSeconds'});
    writetable(summary,fullfile(C.resultsDir, ...
        sprintf('download_summary_%04d.csv',maxStations)));
    save(fullfile(C.resultsDir,sprintf('download_summary_%04d.mat',maxStations)), ...
        'summary','C');
end

function out=download_one(url,destination,station)
    ticValue=tic;
    out=struct('station',station,'status','failed','bytes',0, ...
        'elapsedSeconds',0);
    if isfile(destination)
        info=dir(destination);
        if info.bytes>1000
            out.status='existing'; out.bytes=info.bytes;
            out.elapsedSeconds=toc(ticValue); return;
        end
    end
    for attempt=1:3
        try
            options=weboptions('Timeout',120);
            websave(destination,url,options);
            info=dir(destination);
            out.status='downloaded'; out.bytes=info.bytes;
            out.elapsedSeconds=toc(ticValue); return;
        catch err
            if contains(err.message,'404') || contains(err.message,'Not Found')
                out.status='missing'; break;
            end
            pause(attempt);
        end
    end
    if isfile(destination)
        info=dir(destination);
        if info.bytes<1000, delete(destination); end
    end
    out.elapsedSeconds=toc(ticValue);
end
