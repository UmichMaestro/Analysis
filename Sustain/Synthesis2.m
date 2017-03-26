function [out] = Synthesis2( frq, A, fs, window )
%Synthesis by Jungho!
%   Detailed explanation goes here

if nargin < 4 || isEmpty(window)
    window = 100;
end

samples = fs/window;

harmonics = size(A, 1); % number of harmonics
time = size(A, 2); 
pIncr = 2*pi*frq/fs;

out = zeros(1, time * samples); % 44100 Hz. 100
phase = 0;
for t = 1:(time-1) 
    p = phase;
    for s = 1:samples % each data corresponds to 441 samples.
        p = p + pIncr;
        i = (t-1)*samples+s;
        
        for h = 1:harmonics 
            a0 = A(h,t);
            a1 = A(h,t+1);
            amp = a0 + (a1-a0)*s./samples;    
            out(i) = out(i) + amp * cos(p*h);
        end
    end
    
    phase = p;
end

out = out / rms(out);
