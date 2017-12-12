function [ newA ] = ConvertARModelToAmp( A, start, finish, sustain, model, ar )
    harmonics = size(A,1);
    time = size(A, 2);
    duration = (time - (finish-start)) + sustain;
    newA = zeros(harmonics, duration);
    
    newA(:,1:start-1) = A(:,1:start-1);
    
    if (ar)
        noise = randn(1, duration);
        for i = 1:harmonics 
%             m = model(i,1);
            m = mean([A(i,start),A(i,finish)]);
            g = sqrt(model(i,2));
            a = [1, model(i,3:end)];
            newA(i, start:start+duration-1) = m + filter(a, 1, g*noise);
        end
    else
        for i = 1:harmonics 
%             m = model(i,1);
            m = mean([A(i,start),A(i,finish)]);
            newA(i, start:start+duration-1) = m;
        end
    end

    release = time-finish;
    newA(:,end-release:end) = A(:,end-release:end);
    
end

