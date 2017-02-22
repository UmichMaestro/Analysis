function SplitFile( fileName )
%SPLITFILE Summary of this function goes here
%   Detailed explanation goes here

[segments, fs] = detectVoiced(fileName)
[folder,name] = fileparts(fileName)

for idx = 1:numel(segments)
    s = segments(idx)
    frq = AverageFrequency(s, fs);
    
    outFileName = sprintf('%s/%s-%d-%d.wav',folder,name,idx,frq)
    audiowrite(outFileName,s{1},fs)
end

