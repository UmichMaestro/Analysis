function CrackIowa( url, outputPath ) % http://theremin.music.uiowa.edu/sound%20files/MIS/Woodwinds/oboe/oboe.ff.zip
if nargin < 2 || isempty(outputPath)
    outputPath = '~/Downloads/';
end

[~,n,~] = fileparts(url);

disp('Downloading..........')
output = websave(sprintf('%s%s.zip', outputPath, n), url);

disp('Unzipping............')
path = output(1:end-4);
unzip(output,path)
delete(output)

disp('Analyze!!!!!!!!!!!!!!')
dirResult = dir(path);
for idx = 1:numel(dirResult)
    item = dirResult(idx);
    fileName = item.name;
    [~,~,e] = fileparts(fileName);
    if strcmp(e,'.aiff')
        path = sprintf('%s/%s',item.folder,item.name);
        Process(path)
    end
end

