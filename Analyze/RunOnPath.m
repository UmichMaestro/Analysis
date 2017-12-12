function RunOnPath( path )
%RUNONPATH Summary of this function goes here
%   Detailed explanation goes here
dirResult = dir(path);
for idx = 1:numel(dirResult)
    item = dirResult(idx);
    fileName = item.name;
    [~,~,e] = fileparts(fileName);
    if strcmp(e,'.aiff') || strcmp(e,'.aif')
        path = sprintf('%s/%s',item.folder,item.name);
        Process(path)
    end
end

end

