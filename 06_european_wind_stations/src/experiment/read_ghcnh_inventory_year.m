function inventory = read_ghcnh_inventory_year(filename,targetYear,logFile)
%READ_GHCNH_INVENTORY_YEAR Read monthly record counts for one GHCNh year.

    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open GHCNh inventory: %s',filename); end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fgetl(fid); % header
    capacity = 10000;
    station = strings(capacity,1);
    counts = zeros(capacity,12);
    n = 0; scanned = 0;
    yearText = sprintf('%04d',targetYear);
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        scanned = scanned+1;
        if numel(line)<17 || ~strcmp(line(13:16),yearText), continue; end
        values = sscanf(line(18:end),'%f',12);
        if numel(values)~=12, continue; end
        n=n+1;
        if n>capacity
            capacity=capacity*2;
            station(capacity,1)="";
            counts(capacity,12)=0;
        end
        station(n)=strtrim(string(line(1:11)));
        counts(n,:)=values(:)';
        if mod(n,5000)==0
            ghcnh_log(logFile,'inventory',sprintf( ...
                'found %d station records for %d',n,targetYear));
        end
    end
    station=station(1:n); counts=counts(1:n,:);
    names={'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
    inventory=array2table(counts,'VariableNames',names);
    inventory=addvars(inventory,station,'Before',1,'NewVariableNames','station');
    ghcnh_log(logFile,'inventory',sprintf( ...
        'scanned %d rows and retained %d records for %d',scanned,n,targetYear));
end
