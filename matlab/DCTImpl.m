function y = DCTImpl(x)
    N = length(x);
    if N == 1
        y = x;
    else
        x1 = [x(1:2:N, :); x(N:(-2):1, :)];
        y = FFTImpl(x1, @FFTKernelStandard);
        cosvec = cos(pi*((0:(N-1))')/(2*N));
        sinvec = sin(pi*((0:(N-1))')/(2*N));
        for s2 = 1:size(x, 2)
            y(:, s2) = cosvec.*real(y(:, s2)) + sinvec.*imag(y(:, s2));
        end
        y(1, :) = sqrt(1/N)*y(1, :);
        y(2:N, :) = sqrt(2/N)*y(2:N, :);
    end