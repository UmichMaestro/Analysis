function [out] = Synthesis( frq, A, fs, window )
%Synthesis by Jungho!
%   Detailed explanation goes here

if nargin < 4 || isEmpty(window)
    window = 100;
end

samples = fs/window;

harmonics = size(A, 1); % number of harmonics
time = size(A, 2); 

out = zeros(1, time * samples); % 44100 Hz. 100
for h = 1:harmonics 
    phase = 0;
    for t = 1:(time-1) 
        a0 = A(h,t);
        a1 = A(h,t+1);
        for s = 1:samples % each data corresponds to 441 samples.
            amp = a0 + (a1-a0)*s./samples;
            phase = phase + 2*pi*frq/fs;
            i = (t-1)*samples+s;
            out(i) = out(i) + amp * cos(phase*h);
        end
    end
end

out = out / rms(out);
