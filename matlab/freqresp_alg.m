function freqresp_alg(wave_name, lowpass, dual)
    idwt_kernel = find_kernel(wave_name, 0, dual, 0);
    N = 128;
    n = (0:(N-1))';
    omega = 2*pi*n/N;

    g = zeros(N, 1);
    if lowpass
        g(1) = 1;
    else
        g(2) = 1;
    end
    
    g = idwt_kernel(g, 'per');
    figure();
    plot(omega, abs(fft(g)), 'k-')
end