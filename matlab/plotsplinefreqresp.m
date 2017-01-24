function plotsplinefreqresp(N1, N2)
    N = (N1 + N2)/2;
    h0 = computeQN(N);
    for k=1:(N1/2)
        h0 = conv(h0, [1/4 1/2 1/4]);
    end
    g0 = [1];
    for k=1:(N2/2)
        g0 = conv(g0, [1/4 1/2 1/4]);
    end
    
    L = 100;
    h0 = [h0 zeros(1, L - length(h0))];
    g0 = [g0 zeros(1, L - length(g0))]; 
    omega = 2*pi*(0:(L-1))/L;
    
    subplot(1,2,1);
    plot(omega, abs(fft(h0)))
    axis equal
    
    subplot(1,2,2);
    plot(omega, abs(fft(g0)))
    axis equal