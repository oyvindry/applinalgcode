function freqresp_alg(wave_name, lowpass, dual)
    N = 128;
    n = (0:(N-1))';
    omega = 2*pi*n/N;

    g = zeros(N, 1);
    if lowpass
        g(1) = 1;
    else
        g(2) = 1;
    end
    
    g = idwt_impl(g, wave_name, 1, 'per', 'none', 1, dual, 0, 'time');
    figure();
    plot(omega, abs(fft(g)), 'k-')
end