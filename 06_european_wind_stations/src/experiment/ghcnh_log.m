function ghcnh_log(logFile,stage,message)
%GHCNH_LOG Write the same progress message to the console and a log file.

    folder = fileparts(logFile);
    if ~isfolder(folder), mkdir(folder); end
    stamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    line = sprintf('%s\t%s\t%s\n',stamp,char(stage),char(message));
    fid = fopen(logFile,'a');
    if fid >= 0
        fprintf(fid,'%s',line);
        fclose(fid);
    end
    fprintf('[%s] %s: %s\n',stamp,char(stage),char(message));
end
