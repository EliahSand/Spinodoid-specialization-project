function ensure_dir(folderPath)
%ENSURE_DIR Create folderPath if it does not already exist.
if ~isfolder(folderPath)
    mkdir(folderPath);
end
end
