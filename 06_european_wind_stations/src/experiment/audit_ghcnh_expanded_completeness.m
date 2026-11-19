function [coverage,Y,targetTimes]=audit_ghcnh_expanded_completeness()
%AUDIT_GHCNH_EXPANDED_COMPLETENESS Audit the expanded manifest separately.

    C=ghcnh_trace_config();
    manifest=dir(fullfile(C.dataProcessedDir, ...
        'download_candidates_expanded_*.csv'));
    if isempty(manifest), error('Expanded candidate manifest is absent.'); end
    [~,latest]=max([manifest.datenum]);
    candidates=readtable(fullfile(manifest(latest).folder,manifest(latest).name), ...
        'TextType','string');
    n=height(candidates); targetTimes=make_target_times(C);
    Y=nan(n,numel(targetTimes)); status=strings(n,1);
    nRaw=zeros(n,1); nAccepted=zeros(n,1); nMatched=zeros(n,1);
    ghcnh_log(C.progressLog,'expanded_completeness',sprintf( ...
        'auditing %d files at %d target times',n,numel(targetTimes)));

    for i=1:n
        station=char(candidates.station(i));
        file=fullfile(C.stationDataDir,sprintf('GHCNh_%s_%d.parquet', ...
            station,C.year));
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
            nMatched(i)=sum(isfinite(Y(i,:))); status(i)="audited";
        catch err
            status(i)="parse_failed";
            ghcnh_log(C.progressLog,'expanded_completeness',sprintf( ...
                '%s parse failure: %s',station,err.message));
        end
        if mod(i,25)==0 || i==n
            ghcnh_log(C.progressLog,'expanded_completeness',sprintf( ...
                '%d/%d audited; complete=%d, near-complete=%d',i,n, ...
                sum(nMatched(1:i)==C.expectedTimePoints), ...
                sum(nMatched(1:i)>=C.expectedTimePoints-C.maximumNearCompleteGaps)));
        end
    end

    coverage=candidates; coverage.status=status; coverage.nRaw=nRaw;
    coverage.nAccepted=nAccepted; coverage.nMatched=nMatched;
    coverage.nMissing=C.expectedTimePoints-nMatched;
    coverage.complete=nMatched==C.expectedTimePoints;
    coverage.nearComplete=nMatched>= ...
        C.expectedTimePoints-C.maximumNearCompleteGaps;
    stem=sprintf('expanded_%04d',n);
    writetable(coverage,fullfile(C.dataProcessedDir, ...
        ['completeness_' stem '.csv']));
    save(fullfile(C.dataProcessedDir,['wind_panel_' stem '.mat']), ...
        'C','coverage','Y','targetTimes','-v7.3');
    ghcnh_log(C.progressLog,'expanded_completeness',sprintf( ...
        'finished: %d complete and %d near-complete stations', ...
        sum(coverage.complete),sum(coverage.nearComplete)));
end

function targetTimes=make_target_times(C)
    targetTimes=NaT(C.expectedTimePoints,1,'TimeZone','UTC'); cursor=0;
    for day=0:(C.nDays-1)
        for hour=C.synopticHours
            cursor=cursor+1;
            targetTimes(cursor)=dateshift(C.startTime+days(day), ...
                'start','day')+hours(hour);
        end
    end
end
