function plotboundaryfunctions(N)
    m = 8;
    M = 4*N;

    t = linspace(0, M, M*2^m);
    wave_name = sprintf('db%d', N);

    for k=1:N
        coords = zeros(M*2^m,1);
        coords(k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k);
        plot(t, 2^(m/2)*coords, 'k-');
        title(strcat('\phi^{left}_', sprintf('%u',k-1)))
    end

    for k=1:N
        coords = zeros(M*2^m,1);
        coords(M+k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k+N);
        plot(t, 2^(m/2)*coords, 'k-')
        title(strcat('\psi^{left}_', sprintf('%u',k-1)))
    end
    print('funcsleft','-dpdf')
    cla
    
    for k=1:N
        coords = zeros(M*2^m,1);
        coords(N+k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k);
        plot(t, 2^(m/2)*coords, 'k-')
        title(strcat('\phi_', sprintf('%u',N+k-1)))
    end
    
    for k=1:N
        coords = zeros(M*2^m,1);
        coords(M+N+k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k+N);
        plot(t, 2^(m/2)*coords, 'k-')
        title(strcat('\psi_', sprintf('%u',N+k-1)))
    end
    print('funcsmid','-dpdf')
    cla
    
    for k=1:N
        coords = zeros(M*2^m,1);
        coords(M-N+k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k);
        plot(t, 2^(m/2)*coords, 'k-')
        title(strcat('\phi^{right}_', sprintf('%u',M-N+k-1)))
    end

    for k=1:N
        coords = zeros(M*2^m,1);
        coords(M+M-N+k) = 1;
        coords=IDWTImpl(coords, m, wave_name, 2);
        subplot(2,N,k+N);
        plot(t, 2^(m/2)*coords, 'k-')
        title(strcat('\psi^{right}_', sprintf('%u',M-N+k-1)))
    end
    print('funcsright','-dpdf')
    cla
