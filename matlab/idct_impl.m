function x = idct_impl(y)
    N = size(y, 1);
    if N == 1 
        x = y;
    else
        y(1, :) = y(1, :)/sqrt(1/N);
        y(2:N, :) = y(2:N, :)/sqrt(2/N);
        Q = exp(-pi*1i*((0:(N-1))')/(2*N));
        y1 = zeros(size(y)); y1(1, :) = y(1, :)/Q(1);
        for s2 = 1:size(y, 2)
            y1(2:N, s2) = (y(2:N, s2)-1i*y(N:(-1):2, s2))./Q(2:N);
        end
        y1 = fft_impl(y1, @fft_kernel_standard, 0);
        x = zeros(size(y));
        x(1:2:N, :) = real(y1(1:(N/2), :));
        x(2:2:N, :) = real(y1(N:(-1):(N/2+1), :));
    end
end