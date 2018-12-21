function y = fft_kernel_splitradix(x, forward)
    N = size(x, 1);
    sign = -1;
    if ~forward
        sign = 1;
    end
    if N == 1
        y = x;
    elseif N == 2
        y = [x(1) + x(2); x(1) - x(2)];
    else
        xe  = fft_kernel_splitradix(x(1:(N/2)), forward);
        xo1 = fft_kernel_splitradix(x((N/2 + 1):(3*N/4)), forward);
        xo2 = fft_kernel_splitradix(x((3*N/4 + 1):N), forward);
        G = exp(sign*2*pi*1j*(0:(N/4-1))'/N);
        H = G.*exp(sign*2*pi*1j*(0:(N/4-1))'/(N/2));
        xo1 = G.*xo1;
        xo2 = H.*xo2;
        xo = [xo1 + xo2; -sign*1i*(xo2 - xo1)];
        y = [xe + xo; xe - xo];
    end
end