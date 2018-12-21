function z = mp3_forward_fbt(x)
    N = length(x);
    z = zeros(N,1);
    C = mp3_ctable(); % The analysis window;
    x = flipud(x);
    x = [x; zeros(512-32,1)];
    % The 32x64 matrix M
    M = cos((2*((0:31)')+1)*((0:63)-16)*pi/64);
    
    start = length(x) - 511;
    n = 1;
    for n = 1:(N/32)
        X = x(start:(start + 511));
        Z = C.*X;
        Y = zeros(64, 1);
        for j = 0:7
            Y = Y + Z((64*j + 1):(64*(j + 1)));
        end
        z((1+(n-1)*32):(n*32)) = M*Y;
        start = start - 32;
    end















