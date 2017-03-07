function x = DFTImpl8(x)
    N = size(x, 1);
    for n = 1:8:N
        x(n:(n+7), :) = fft(x(n:(n+7), :));
    end
end