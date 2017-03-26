start = 1;
finish = 5;
loop = 10;

fs = 44100;
window = 100;

%[frq, A] = OpenBinary('Horn.ff.C4B4-1-261.msm');
[frq, A] = OpenBinary('/Users/bangtoven/Google Drive/2017 Winter/MAESTRO/WINTER 2017 Files/Synth/_Maestro Synth Model/woodwind/oboe/oboe.ff.C4B4-2-277.msm');
samples = fs/window;

harmonics = size(A, 1); % number of harmonics
time = size(A, 2);
pIncr = 2*pi*frq/fs;

out = zeros(1, (time + (loop-1)*(finish-start))*samples); % 44100 Hz. 100

phase = 0;
outIndex = 1;

[out, phase, outIndex] = SynthInRegion(1, start, phase, pIncr, samples, harmonics, outIndex, A, out);
for l = 1:loop
    [out, phase, outIndex] = SynthInRegion(start, finish, phase, pIncr, samples, harmonics, outIndex, A, out);
end
[out, ~, ~] = SynthInRegion(finish, time, phase, pIncr, samples, harmonics, outIndex, A, out);
out = out / rms(out) / 32;
player = audioplayer(out, fs);
play(player)

function [out, phase, outIndex] = SynthInRegion(start, finish, phase, pIncr, samples, harmonics, outIndex, A, out)
    for t = start:finish-1
        p = phase;
        for s = 1:samples % each data corresponds to 441 samples.
            p = p + pIncr;
            i = (outIndex-1)*samples + s;

            for h = 1:harmonics
                a0 = A(h,t);
                a1 = A(h,t+1);
                amp = a0 + (a1-a0)*s./samples;
                out(i) = out(i) + amp * cos(p*h);
            end
        end

        phase = p;
        outIndex = outIndex + 1;
    end
end
