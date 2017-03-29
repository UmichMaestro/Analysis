start = 30;
finish = 30;
loop = 5;


fs = 44100;
window = 100;

%[frq, A] = OpenBinary('Horn.ff.C4B4-1-261.msm');
[frq, A] = OpenBinary('/Users/bangtoven/Google Drive/2017 Winter/MAESTRO/WINTER 2017 Files/Synth/_Maestro Synth Model/woodwind/oboe/oboe.ff.C4B4-2-277.msm');

samples = fs/window;

harmonics = size(A, 1); % number of harmonics
time = size(A, 2); 

out = zeros(1, time*samples + (loop-1)*(time-start-finish)*samples + 1000000); % 44100 Hz. 100
for h = 1:harmonics 
    phase = 0;
    
    [out, phase] = Calculate(A, frq, fs, samples, h, phase, out, 1, 1+start-1, 1);
    
    startIndex = start-1;
    for l = 1:loop
        [out, phase] = Calculate(A, frq, fs, samples, h, phase, out, 1+start, time-1-finish-1, startIndex);
        startIndex = startIndex + (time-start-finish-3);
    end
    
    [out, phase] = Calculate(A, frq, fs, samples, h, phase, out, time-1-finish, time-1, startIndex);
end

out = out / rms(out);
out = out / 32;
player = audioplayer(out, fs);
playblocking(player)


function [out, phase] = Calculate(A, frq, fs, samples, h, phase, out, start, finish, startIndex)
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