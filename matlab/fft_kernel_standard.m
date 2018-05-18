% Compute the DFT of the column vector x using the FFT algorithm
function y = fft_kernel_standard(x, forward)
    N = size(x, 1);
    sign = -1;
    if ~forward
        sign = 1;
    end
    if N == 1 
        y = x;
    else
        xe = fft_kernel_standard(x(1:(N/2)), forward);
        xo = fft_kernel_standard(x((N/2+1):N), forward);
        D = exp(sign*2*pi*1j*(0:(N/2-1))'/N);
        xo = xo.*D;
        y = [ xe + xo; xe - xo];
    end
end