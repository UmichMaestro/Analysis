function [freq, A] = Analysis( signal, fs )
    [frq_path, ac_path, t, amp] = Praat(signal{1}, fs);

    upSynthFreq = 10000;
    iGood = find(frq_path > 0);
    mnF0 = mean(frq_path(iGood));
    nHarm = round(upSynthFreq/mnF0);    % adjust nHarm to reflect f0 - synthesize up to 8...10kHz?

    freq = round(mnF0)
    
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
    for iSmp = 1:length(t)
        tmp = round(t(iSmp))+ctrIdx;
        iLoc = find(tmp >= 1 & tmp <= t(end));
        idx = tmp(iLoc);
        for iHarm = 1:nHarm
            A(iHarm,iSmp) = abs(dot(exp(j*2*pi*iHarm*frq_path(iSmp)*idx/fs),signal{1}(idx)'));
        end
    end