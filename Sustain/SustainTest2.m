start = 0;
finish = 0;
loop = 1;


fs = 44100;
window = 100;

[frq, A] = OpenBinary('Horn.ff.C4B4-1-261.msm');

samples = fs/window;

harmonics = size(A, 1); % number of harmonics
time = size(A, 2); 

out = zeros(1, time * samples + (loop-1)*(time-start-finish)*samples+10000000); % 44100 Hz. 100
    
[out, phase] = Calculate(A, frq, fs, samples, h, harmonics, phase, out, 1, 1+start-1, 1);

startIndex = 1+start-1;
for l = 1:loop
    [out, phase] = Calculate(A, frq, fs, samples, h, harmonics, phase, out, 1+start, time-1-finish-1, startIndex);
    startIndex = startIndex + (time-start-finish-1);
end

[out, phase] = Calculate(A, frq, fs, samples, h, harmonics, phase, out, time-1-finish, time-1, startIndex);

out = out / rms(out);
out = out / 32;
player = audioplayer(out, fs);
play(player)


function [out, phase] = Calculate(A, frq, fs, samples, h, harmonics, phase, out, start, finish, startIndex)
for h = 1:harmonics
    phase = 0;
    for t = start:finish
        
        a0 = A(h,t);
        a1 = A(h,t+1);
        for s = 1:samples % each data corresponds to 441 samples.
            amp = a0 + (a1-a0)*s./samples;
            phase = phase + 2*pi*frq/fs;
            i = ((t-1)+startIndex)*samples+s; % check this later.
            out(i) = out(i) + amp * cos(phase*h);
        end
    end
    
end
end