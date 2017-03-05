function Spectrogram( sig, fs, outFileName )
%SPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here
freqRes = 20;       % frequency resolution [Hz] (determines duration of temporal window)
nFFTFac = 2;        % size of FFT relative to size of temporal window
ovlpFac = 0.8;      % proportion of overlap between consecutive temporal windows
beaut.range = 100;        % range of intensities to display [dB]
beaut.fMin = 0;           % minimum frequency to display [Hz]
beaut.fMax = 10000;
nWin = floor(fs/freqRes);           % size of temporal window
nfft = 2^ceil(log2(nFFTFac*nWin));  % FFT size
ovlp = floor(ovlpFac*nWin);
win = hamming(nWin);

[y,f,t,p] = spectrogram(sig(:,1),win,ovlp,nfft,fs,'yaxis');

% Display spectrum in decibels over desired ranges
pdB = 10*log10(abs(p));    % Transform time-frequency spectrum (power) into decibels

vMax = max(max(pdB));      % Threshold pdB values to fit within intensity range
iLoc = find(pdB <= vMax - beaut.range);
pdB(iLoc) = vMax - beaut.range;

idf = find(f <= beaut.fMax & f >= beaut.fMin);  % Determine range of frequencies to display

% ax1 = axes;
surf(t,f(idf),pdB(idf,:),'EdgeColor','none');
% view(ax1,0,90);
% axis(ax1,'xy'); 
% axis(ax1,'tight'); 
view(2);
print(strcat(outFileName,'.png'),'-dpng')
% cla(ax1);

end

