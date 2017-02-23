function Process( fileName )
%SPLITFILE Summary of this function goes here
%   Detailed explanation goes here

[segments, fs] = SegmentSound(fileName);
[folder,name] = fileparts(fileName);

for idx = 1:numel(segments)
    s = segments(idx);
    [frq, A] = Analysis(s, fs);
    
    outFileName = sprintf('%s/%s-%d-%d',folder,name,idx,round(frq));
    audiowrite(strcat(outFileName,'.wav'),s{1},fs)
    csvwrite(strcat(outFileName,'.csv'), A')
    plot(A')
    print(strcat(outFileName,'.png'),'-dpng')
    
    % synthesize
    outsig = Synthesis(frq, A, fs);
    % normalize
    outsig = rms(s{1}) / rms(outsig) * outsig;
    audiowrite(sprintf('%s-synth.wav', outFileName), outsig, fs)
end
    

        