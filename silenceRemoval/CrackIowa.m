function CrackIowa( url )
%CRACKIOWA Summary of this function goes here
%   Detailed explanation goes here

% http://theremin.music.uiowa.edu/sound%20files/MIS/Woodwinds/oboe/oboe.ff.zip
[~,n,~] = fileparts(url);

disp('Downloading..........')
output = websave(sprintf('~/Downloads/%s.zip', n), url);

disp('Unzipping............')
path = output(1:end-4);
unzip(output,path)
delete(output)

disp('Analyze!!!!!!!!!!!!!!')
SplitFilesInFolder(path)

end

