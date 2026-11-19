function [coverage,Y,targetTimes] = audit_ghcnh_completeness(maxStations)
%AUDIT_GHCNH_COMPLETENESS Match quality-controlled wind to a regular grid.

    C=ghcnh_trace_config();
    if nargin<1 || isempty(maxStations), maxStations=C.downloadCandidateCount; end
    candidates=readtable(fullfile(C.dataProcessedDir,'download_candidates.csv'), ...
        'TextType','string');
    maxStations=min(maxStations,height(candidates));
    candidates=candidates(1:maxStations,:);
    targetTimes=make_target_times(C);
    Y=nan(maxStations,numel(targetTimes));
    status=strings(maxStations,1);
    nRaw=zeros(maxStations,1); nAccepted=zeros(maxStations,1);
    nMatched=zeros(maxStations,1);

    ghcnh_log(C.progressLog,'completeness',sprintf( ...
        'auditing %d files at %d target times',maxStations,numel(targetTimes)));
    for i=1:maxStations
        station=char(candidates.station(i));
        file=fullfile(C.stationDataDir,sprintf('GHCNh_%s_%d.parquet',station,C.year));
        if ~isfile(file), status(i)="file_missing"; continue; end
        try
            tbl=parquetread(file,'SelectedVariableNames', ...
                {'DATE','wind_speed','wind_speed_Quality_Code'});
            nRaw(i)=height(tbl);
            time=datetime(string(tbl.DATE),'InputFormat', ...
                "yyyy-MM-dd'T'HH:mm:ss",'TimeZone','UTC');
            speed=double(tbl.wind_speed);
            quality=string(tbl.wind_speed_Quality_Code);
            accepted=isfinite(speed) & quality=="1" & ...
                time>=targetTimes(1)-minutes(C.timeToleranceMinutes) & ...
                time<=targetTimes(end)+minutes(C.timeToleranceMinutes);
            time=time(accepted); speed=speed(accepted);
            nAccepted(i)=numel(speed);
            for t=1:numel(targetTimes)
                [distance,index]=min(abs(minutes(time-targetTimes(t))));
                if ~isempty(distance) && distance<=C.timeToleranceMinutes
                    Y(i,t)=speed(index);
                end
            end
            nMatched(i)=sum(isfinite(Y(i,:)));
            status(i)="audited";
        catch err
            status(i)="parse_failed";
            ghcnh_log(C.progressLog,'completeness',sprintf( ...
                '%s parse failure: %s',station,err.message));
        end
        if mod(i,C.progressEveryStations)==0 || i==maxStations
            ghcnh_log(C.progressLog,'completeness',sprintf( ...
                '%d/%d audited; complete=%d, near-complete=%d',i,maxStations, ...
                sum(nMatched(1:i)==C.expectedTimePoints), ...
                sum(nMatched(1:i)>=C.expectedTimePoints-C.maximumNearCompleteGaps)));
        end
    end
    coverage=candidates;
    coverage.status=status; coverage.nRaw=nRaw; coverage.nAccepted=nAccepted;
    coverage.nMatched=nMatched;
    coverage.nMissing=C.expectedTimePoints-nMatched;
    coverage.complete=nMatched==C.expectedTimePoints;
    coverage.nearComplete=nMatched>=C.expectedTimePoints-C.maximumNearCompleteGaps;
    if ~isfolder(C.dataProcessedDir), mkdir(C.dataProcessedDir); end
    writetable(coverage,fullfile(C.dataProcessedDir, ...
        sprintf('completeness_%04d.csv',maxStations)));
    save(fullfile(C.dataProcessedDir,sprintf('wind_panel_%04d.mat',maxStations)), ...
        'C','coverage','Y','targetTimes','-v7.3');
    ghcnh_log(C.progressLog,'completeness',sprintf( ...
        'finished: %d complete and %d near-complete stations', ...
        sum(coverage.complete),sum(coverage.nearComplete)));
end

function targetTimes=make_target_times(C)
    targetTimes=NaT(C.expectedTimePoints,1,'TimeZone','UTC');
    cursor=0;
    for day=0:(C.nDays-1)
        for hour=C.synopticHours
            cursor=cursor+1;
            targetTimes(cursor)=dateshift(C.startTime+days(day),'start','day')+hours(hour);
        end
    end
end
