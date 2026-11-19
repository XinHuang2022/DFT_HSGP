function summary=download_ghcnh_expanded_candidates()
%DOWNLOAD_GHCNH_EXPANDED_CANDIDATES Download remaining expanded-manifest files.

    C=ghcnh_trace_config();
    manifest=dir(fullfile(C.dataProcessedDir, ...
        'download_candidates_expanded_*.csv'));
    if isempty(manifest)
        prepare_ghcnh_expanded_manifest();
        manifest=dir(fullfile(C.dataProcessedDir, ...
            'download_candidates_expanded_*.csv'));
    end
    [~,latest]=max([manifest.datenum]);
    candidates=readtable(fullfile(manifest(latest).folder,manifest(latest).name), ...
        'TextType','string');
    n=height(candidates);
    if ~isfolder(C.stationDataDir), mkdir(C.stationDataDir); end
    ghcnh_log(C.progressLog,'expanded_download',sprintf( ...
        'starting/resuming %d station files with %d workers',n,C.downloadWorkers));

    pool=gcp('nocreate');
    if isempty(pool) || pool.NumWorkers~=C.downloadWorkers
        if ~isempty(pool), delete(pool); end
        pool=parpool('local',C.downloadWorkers);
    end
    futures=parallel.FevalFuture.empty(n,0);
    for i=1:n
        station=char(candidates.station(i));
        url=sprintf(C.stationParquetUrl,C.year,station,C.year);
        destination=fullfile(C.stationDataDir, ...
            sprintf('GHCNh_%s_%d.parquet',station,C.year));
        futures(i)=parfeval(pool,@download_one,1,url,destination,station);
    end

    station=strings(n,1); status=strings(n,1);
    bytes=zeros(n,1); elapsed=zeros(n,1);
    for completed=1:n
        [index,out]=fetchNext(futures);
        station(index)=string(out.station); status(index)=string(out.status);
        bytes(index)=out.bytes; elapsed(index)=out.elapsedSeconds;
        if mod(completed,25)==0 || completed==n
            ghcnh_log(C.progressLog,'expanded_download',sprintf( ...
                '%d/%d complete; ok=%d, missing=%d, failed=%d; %.1f MB', ...
                completed,n,sum(status=="downloaded" | status=="existing"), ...
                sum(status=="missing"),sum(status=="failed"),sum(bytes)/2^20));
        end
    end
    summary=table(station,status,bytes,elapsed, ...
        'VariableNames',{'station','status','bytes','elapsedSeconds'});
    writetable(summary,fullfile(C.resultsDir,sprintf( ...
        'download_summary_expanded_%04d.csv',n)));
    save(fullfile(C.resultsDir,sprintf( ...
        'download_summary_expanded_%04d.mat',n)),'summary','C');
end

function out=download_one(url,destination,station)
    ticValue=tic;
    out=struct('station',station,'status','failed','bytes',0, ...
        'elapsedSeconds',0);
    if isfile(destination)
        info=dir(destination);
        if info.bytes>1000
            out.status='existing'; out.bytes=info.bytes;
            out.elapsedSeconds=toc(ticValue); return
        end
    end
    for attempt=1:3
        try
            websave(destination,url,weboptions('Timeout',120));
            info=dir(destination); out.status='downloaded';
            out.bytes=info.bytes; out.elapsedSeconds=toc(ticValue); return
        catch err
            if contains(err.message,'404') || contains(err.message,'Not Found')
                out.status='missing'; break
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
