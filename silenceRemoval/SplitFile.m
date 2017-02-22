function SplitFile( fileName )
%SPLITFILE Summary of this function goes here
%   Detailed explanation goes here

[segments, fs] = detectVoiced(fileName)
[folder,name] = fileparts(fileName)

for idx = 1:numel(segments)
    s = segments(idx)
    [frq, A] = Analysis(s, fs);
    
    outFileName = sprintf('%s/%s-%d-%d',folder,name,idx,frq)
    audiowrite(strcat(outFileName,'.wav'),s{1},fs)
    csvwrite(strcat(outFileName,'.csv'), A')
    plot(A')
    print(strcat(outFileName,'.png'),'-dpng')
end
    

        