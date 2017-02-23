function [freq, A] = Analysis( signal, fs )
    [frq_path, ac_path, time, amp] = Praat(signal{1}, fs);

    upSynthFreq = 10000;
    iGood = find(frq_path > 0);
    mnF0 = mean(frq_path(iGood));
    nHarm = round(upSynthFreq/mnF0);    % adjust nHarm to reflect f0 - synthesize up to 8...10kHz?

    freq = mnF0
    
    if mnF0 > 220
        frameSz = 2048;   
    elseif mnF0 > 110
        frameSz = 2*2048;
    elseif mnF0 > 55
        frameSz = 4*2048;
    else
        frameSz = 8*2048;
    end

    ctrIdx = -frameSz/2+1:frameSz/2;
    pitchOff = 1;%//;2^(11/12);
    for iSmp = 1:length(time)
        tmp = round(time(iSmp))+ctrIdx;
        iLoc = find(tmp >= 1 & tmp <= time(end));
        idx = tmp(iLoc);
        for iHarm = 1:nHarm
            A(iHarm,iSmp) = abs(dot(exp(j*2*pi*iHarm*frq_path(iSmp)*idx/fs),signal{1}(idx)'));
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [frq_path, ac_path, t, amp] = Praat(in, fs, disp_progress, ...
    time_step, min_frq, max_frq, VoiceThresh, SilThresh, OctaveCost,...
    VUnvCost, OctJumpCost, Hypotheses)

%function [frq_path, ac_path, t, amp] = praat_pd_mod(in, fs, disp_progress, ...
%    time_step, min_frq, max_frq, VoiceThresh, SilThresh, OctaveCost,...
%    VUnvCost, OctJumpCost, Hypotheses)
%
% PRAAT_PD_MOD  -- An autocorrelation-based pitch tracker.
%   Note: Parameters can be set to their default value by setting them to
%         empty ([]) or by truncating the parameter list.
%
% Input Parameters:
%   in            -- Single-channel input waveform
%   fs            -- Sampling frequency (samples per second)
%   disp_progress -- 0 or 1; Use a textual waitbar for progress? (Default: 1)
%   time_step     -- EITHER the (scalar) time step in seconds (Default: 0.01)
%                      OR a vector of samples about which to center frames
%   min_frq       -- Minimum frequency in Hertz (Default: 65 Hz)
%                      Should be as high as possible to speed calculation
%   max_frq       -- Maximum frequency in Hertz (Default: 1000 Hz)
%   VoiceThresh   -- Roughly determines the minimum correlation that will be
%                      labeled as "voiced" (Default: 0.5)
%   SilThresh     -- Roughly determines which frames will be labeled as silent
%                      (Default 0.02)
%   OctaveCost    -- Weighting factor to encourage selection of low-lag peaks
%                      (Default 0.02);
%   VUnvCost      -- Cost of a voiced-unvoiced transition (Default: 0.2)
%   OctJumpCost   -- Cost of an octave jump (Default: 0.2)
%   Hypotheses    -- Number of candidate peaks to select at each frame
%                      (Default: 6)
% Output Parameters:
%   frq_pth       -- Output fundamental frequency samples in Hertz
%   ac_pth        -- Autocorrelation samples, effectively a confidence value
%   t             -- Vector of samples about which frames are centered
%   amp           -- The windowed RMS amplitude of output frames

% Based on:
%   P. Boersma, "Accurate short-term analysis of the fundamental frequency
%      and the harmonics-to-noise ratio of a sampled sound," in Proceedings 
%      of the Institute of Phonetic Sciences of the University of Amsterdam, 
%      vol. 17, 1993, pp. 97--110.
%
% With some modifications as described in: 
%   M. A. Bartsch, "Automatic assesment of the spasmodic voice," Communications
%      and Signal Processing Laboratory (CSPL) Tech. Report #339, EECS 
%      Department, University of Michigan, 2003.

%   Author: M. Bartsch (mbartsch@umich.edu)
%   $Revision: 1.00 $  $Date: 2004/04/08 16:30:00 $


if nargin < 2 | isempty(fs)
    fs = 44100;
end
if nargin < 3 | isempty(disp_progress)
    disp_progress = 0;
end
if nargin < 4 | isempty(time_step)
    time_step = 0.01;
end
if nargin < 5 | isempty(min_frq)
    min_frq = 65;
end
if nargin < 6 | isempty(max_frq)
    max_frq = 1000;
end
if nargin < 7 | isempty(VoiceThresh)
%    VoiceThresh = 0.4;
    VoiceThresh = 0.5;
end
if nargin < 8 | isempty(SilThresh)
%    SilThresh = 0.05;
    SilThresh = 0.02;
end
if nargin < 9 | isempty(OctaveCost)
    OctaveCost = 0.02;
end
if nargin < 10 | isempty(VUnvCost)
    VUnvCost = 0.2;
end
if nargin < 11 | isempty(OctJumpCost)
    OctJumpCost = 0.2;
end
if nargin < 12 | isempty(Hypotheses)
    Hypotheses = 6;
end

ac_mat = [];

in = in(:);

%z = fft(in);
%ind = round(length(z)*.475):ceil(length(z)/2);
%z(ind) = z(ind).*linspace(1,0,length(ind))';
%z(length(z)-ind) = z(length(z)-ind).*linspace(1,0,length(ind))';
%z = [z(1:floor(length(z)/2)); 0; z(ceil(length(z)/2):length(z))];
%in = ifft(z);
%
%  The above is a "preprocessing step" given by Boersma... 
%     I'm not sure that it's worth the computation time.  MB

%  In case we want to use longer windows than would be given by the minimum
%   frequency, we can specify a minimum frequency exclusively for the purposes
%   of "window calculation" by setting min_frq to a vector as:
%           min_frq = [act_min_frq, wnd_min_frq]
if length(min_frq) > 1
    wnd_min_frq = min_frq(2);
    min_frq = min_frq(1);
else
    wnd_min_frq = min_frq;
end

window_size = min([round(3*fs/wnd_min_frq), length(in)]);
fft_size = 2^(ceil(log2(window_size*1.5)));
fft_pad = fft_size - window_size;

global_peak = max(abs(in));

window = hanning(window_size);
window_ac = fft([window; zeros(fft_pad,1)]);
window_ac = window_ac.*conj(window_ac);
window_ac = real(ifft(window_ac));
window_ac = window_ac./window_ac(1);

wnd_fft = conj(fft(window,fft_size));

%  Allow the user to specify samples about which to center data windows
if length(time_step) > 1
    t = time_step(:) - floor(window_size/2);
    t(t < 1) = 1;
    t(t > length(in)-window_size+1) = length(in)-window_size+1;
else  
    samp_step = round(time_step*fs);
    t = (1:samp_step:length(in)-window_size)';
end

amp = zeros(length(t),1);
lags = zeros(length(t),Hypotheses);
ac_vals = zeros(length(t),Hypotheses-1);
scores = zeros(length(t),Hypotheses);

wb = waitbar_text(0,'Computing pitch-track.',~disp_progress);
counter = 1;
for time = t'
    frame = in(time:time+window_size-1);
    amp(counter) = sqrt(mean((window.*(frame-mean(frame))).^2));
    
    spec = fft((frame - mean(frame)).*window,fft_size);
    spec2 = fft((frame - mean(frame)).^2.*window,fft_size);

    num = real(ifft(spec.*conj(spec)));
    den1 = real(ifft(spec2.*wnd_fft));
    den2 = den1([1 end:-1:2]);
    den = sqrt(den1.*den2);
    
    ac = repmat(inf,size(num));
    nz = den ~= 0;
    ac(nz) = num(nz)./den(nz);
       
    [lags(counter,1:Hypotheses-1),ac_vals(counter,1:Hypotheses-1)] = ...
        peak_search(ac,Hypotheses-1,fs/max_frq,fs/min_frq,OctaveCost);

    loc_peak = max(frame);
    scores(counter,Hypotheses) = VoiceThresh + ...
        max([0, 2-(loc_peak/global_peak)/(SilThresh/(1+VoiceThresh))]);
    scores(counter,1:Hypotheses-1) = ac_vals(counter,1:Hypotheses-1) ...
        - OctaveCost*log2(min_frq*lags(counter,1:Hypotheses-1)/fs);
    
    wb = waitbar_text(time/t(end),wb,~disp_progress);
    counter = counter + 1;
end
waitbar_text([],[],~disp_progress);
scores(isinf(scores)) = -inf;

frqs = repmat(0,size(lags));
nz = lags ~= 0;
frqs(nz) = fs./lags(nz);


% Implement Viterbi path finding

last_cost = -scores(1,:);
this_cost = zeros(1,size(frqs,2));
paths = zeros(size(frqs));
for time = 2:size(frqs,1)
    for this = 1:size(frqs,2)
        temp_cost = zeros(1,size(frqs,2));
        for that = 1:size(frqs,2)
            if xor(frqs(time-1,that) == 0,frqs(time,this) == 0)
                temp_cost(that) = last_cost(that) + VUnvCost - scores(time,this);
            elseif and(frqs(time-1,that),frqs(time,this))
                cst = OctJumpCost*abs(log2(frqs(time-1,that)/frqs(time,this)));
                temp_cost(that) = last_cost(that) + cst - scores(time,this);
            else
                temp_cost(that) = last_cost(that) - scores(time,this);
            end
        end
        [this_cost(this),paths(time,this)] = min(temp_cost);
    end
    last_cost = this_cost;
    this_cost = zeros(size(this_cost));
end

frq_path = zeros(size(frqs,1),1);
ac_path = zeros(size(frq_path));
[mn,last_choice] = min(last_cost);

for i=length(frq_path):-1:1
    if last_choice == Hypotheses
        ac_path(i) = 0;
    else
        ac_path(i) = ac_vals(i,last_choice);
    end
    frq_path(i) = frqs(i,last_choice);
    last_choice = paths(i,last_choice);
end    

t = t + floor(window_size/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handle = waitbar_text(frac,handle,no_waitbar)
%
% A text-based replacement for MATLAB's waitbar.
% Won't do anything if no_waitbar is equal to 1.

if ~isempty(no_waitbar) & ~no_waitbar
    if nargin < 1
        fprintf(2,'\n');
        return
    end
    L = 60;
    if ischar(handle)
        fprintf(2,'%s\n',handle);
        fprintf(2,'|------------25%%------------50%%------------75%%-------------|\n');
        tmp = 0;
        while tmp < frac
            fprintf(2,'*');
            tmp = tmp + 1/L;
        end
        handle = frac;
    elseif frac > handle
        while handle < frac
            fprintf(2,'*');
            handle = handle + 1/L;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peaks,ampl] = peak_search(x, num_peaks, ind_min, ind_max, OctCost)
%
% Iteratively finds up to 'num_peaks' peaks in the signal 'x'.  
% - Uses parabolic interpolation to refine the location and height of the peak.
% - The search range is restricted by 'ind_min' and 'ind_max'.
% - OctCost is equivalent to Boersma's OctaveCost, and defines to what extent
%     the algorithm prefers peaks with low lag (high frequency).

peak_sep = ind_min;

% Need (ind_min+1) and (ind_max+1) because ind_min, ind_max assume a that
%   DC is the 0th element in the DFT array
ind_min = round(ind_min+1);
ind_max = ind_max+1;

peaks = ones(1,num_peaks);
ampl = ones(1,num_peaks)*-inf;

for fst_neg=1:ind_max
    if x(fst_neg) < 0
        break
    end
end
ind_min = max([fst_neg ind_min]);

inds = find(diff(sign(diff(x(round(ind_min):round(ind_max))))) == -2) + ind_min;

cost = -(x(inds) - OctCost*log2(inds/ind_max));
[srt,ind] = sort(cost);
ind = inds(ind);
ind(x(ind) < 0) = [];
ind(isinf(x(ind))) = [];

for ct = 1:num_peaks
    if isempty(ind)
        return
    end
    
    xval = [ind(1)-2; ind(1)-1; ind(1)];
    A = [xval.^2 xval [1; 1; 1]];
    b = x(xval+1);
    y = A\b;
    peaks(ct) = -y(2)/2/y(1);
    ampl(ct) = y(1)*(peaks(ct))^2 + y(2)*peaks(ct) + y(3);
    
    ind(abs(ind(1) - ind) < peak_sep) = [];
end

