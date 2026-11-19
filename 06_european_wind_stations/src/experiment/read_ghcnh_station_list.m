function stations = read_ghcnh_station_list(filename)
%READ_GHCNH_STATION_LIST Parse the official fixed-width GHCNh station list.

    arguments
        filename (1,:) char
    end
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open station list: %s',filename); end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    raw = textscan(fid,'%s','Delimiter','\n','Whitespace','');
    raw = string(raw{1});
    n = numel(raw);

    station = strings(n,1); latitude = nan(n,1); longitude = nan(n,1);
    elevation = nan(n,1); state = strings(n,1); name = strings(n,1);
    wmo = strings(n,1); country = strings(n,1); network = strings(n,1);
    keep = false(n,1);
    for i = 1:n
        line = char(raw(i));
        if numel(line) < 81, continue; end
        station(i) = strtrim(string(line(1:11)));
        latitude(i) = str2double(line(13:20));
        longitude(i) = str2double(line(22:30));
        elevation(i) = str2double(line(32:37));
        if numel(line)>=40, state(i)=strtrim(string(line(39:40))); end
        if numel(line)>=71, name(i)=strtrim(string(line(42:71))); end
        if numel(line)>=85, wmo(i)=strtrim(string(line(81:85))); end
        if numel(line)>=88, country(i)=strtrim(string(line(87:88))); end
        network(i)="";
        keep(i)=station(i)~="" && isfinite(latitude(i)) && isfinite(longitude(i));
    end
    stations = table(station(keep),latitude(keep),longitude(keep), ...
        elevation(keep),state(keep),name(keep),wmo(keep),country(keep), ...
        network(keep),'VariableNames',{'station','latitude','longitude', ...
        'elevation','state','name','wmo','country','network'});
end
