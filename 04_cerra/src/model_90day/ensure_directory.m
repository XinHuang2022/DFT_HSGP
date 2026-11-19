function ensure_directory(pathName)
%ENSURE_DIRECTORY Create a directory if it does not already exist.

    if ~isfolder(pathName)
        mkdir(pathName);
    end
end


