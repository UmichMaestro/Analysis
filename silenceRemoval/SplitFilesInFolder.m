function SplitFilesInFolder( directory, extension )
%SPLITFILESINFOLDER Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2 || isempty(extension)
    extension = '.aiff'
end

dirResult = dir(directory)

for idx = 1:numel(dirResult)
    item = dirResult(idx)
    fileName = item.name
    [~,~,e] = fileparts(fileName)
    if strcmp(e,extension)
        path = sprintf('%s/%s',item.folder,item.name)
        SplitFile(path)
    end
end
