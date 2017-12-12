% urls = textread('iowa_urls.txt', '%s');
% outPath = '~/Downloads/';
% for idx = 1:length(urls)
%     disp(idx)
%     url = urls{idx};
%     CrackIowa(url, outPath);
% end

function CrackIowa( url, outputPath ) % http://theremin.music.uiowa.edu/sound%20files/MIS/Woodwinds/oboe/oboe.ff.zip
if nargin < 2 || isempty(outputPath)
    outputPath = '~/Downloads/';
end

[~,n,~] = fileparts(url);

disp('Downloading..........')
zipfile = websave(sprintf('%s%s.zip', outputPath, n), url);

disp('Unzipping............')
path = zipfile(1:end-4);
unzip(zipfile,path)
% delete(zipfile)

disp('Analyze!!!!!!!!!!!!!!')
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

function Process( fileName )
[segments, fs] = SegmentSound(fileName);
[folder,name] = fileparts(fileName);

for idx = 1:numel(segments)
    s = segments(idx);
    [frq, A] = Analysis(s, fs);
    
    if isnan(frq)
        outFileName = sprintf('%s/%s-%d',folder,name,idx);
        audiowrite(sprintf('%s-ERR_orig.wav', outFileName),s{1},fs)
        continue
    end
    
    % synthesize & normalize
    outsig = Synthesis(frq, A, fs);
    outsig = rms(s{1}) / rms(outsig) * outsig;
    
    % output
    outFileName = sprintf('%s/%s-%d-%d',folder,name,idx,round(frq));
    % 0. data
    csvwrite(strcat(outFileName,'.csv'), A')
    % 1. audio
    audiowrite(sprintf('%s-wav_orig.wav', outFileName),s{1},fs)
    audiowrite(sprintf('%s-wav_synth.wav', outFileName), outsig, fs)
    % 2. graphs
    plot(A')
    print(sprintf('%s-amp_orig.png', outFileName),'-dpng')
    Spectrogram(s{1}, fs, strcat(outFileName,'-sp_orig'));
    Spectrogram(outsig', fs, strcat(outFileName,'-sp_synth'));
    
end
