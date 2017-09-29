function [ output ] = GenerateModel( A, start, finish, order )
    num = size(A, 1);
    output = zeros(num, order+2);
    for i = 1:num
        sustain = A(i,start:finish); % data of i-th harmonics in sustain portion
        m = mean(sustain);
        r = detrend(sustain, 'linear');
        [a,g] = lpc(r, order);
        output(i, :) = [m, g, a(2:end)];
    end
end
