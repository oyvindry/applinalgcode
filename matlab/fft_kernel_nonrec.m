function y = fft_kernel_nonrec(x, forward)
    N = size(x, 1);
    sign = -1;
    if ~forward
        sign = 1;
    end
    D = exp(sign*2*pi*1i*(0:(N/2 - 1))'/N);
    nextN = 1;
    while nextN < N
        k = 1;
        while k <= N
            xe = x(k:(k + nextN - 1));
            xo = x((k + nextN):(k + 2*nextN - 1)); 
            xo = xo.*D(1:(N/(2*nextN)):(N/2));
            x(k:(k + 2*nextN - 1)) = [xe + xo; xe - xo];
            k = k + 2*nextN;
        end
        nextN = nextN*2;
    end
    y = x;
end