function plotsymlets(N);
    m = 8;
    M = 2*N;

    t = linspace(0, M, M*2^m);
    wave_name = sprintf('sym%d', N);
    
    k = 1;
    coords = zeros(M*2^m,1);
    coords(N) = 1;
    coords=IDWTImpl(coords, m, wave_name, 'per');
    subplot(2,1,k);
    plot(t, 2^(m/2)*coords, 'k-')
    title(strcat('\phi_', sprintf('%u',N+k-1)))
    
    coords = zeros(M*2^m,1);
    coords(M+N) = -1;
    coords=IDWTImpl(coords, m, wave_name, 'per');
    subplot(2,1,k+1);
    plot(t, 2^(m/2)*coords, 'k-')
    title(strcat('\psi_', sprintf('%u',N+k-1)))

end
