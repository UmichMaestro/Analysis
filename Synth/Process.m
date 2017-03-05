function Process( fileName )
%SPLITFILE Summary of this function goes here
%   Detailed explanation goes here

[segments, fs] = SegmentSound(fileName);
[folder,name] = fileparts(fileName);

for idx = 1:numel(segments)
    s = segments(idx);
    [frq, A] = Analysis(s, fs);
    
    if isnan(frq)
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
    

        