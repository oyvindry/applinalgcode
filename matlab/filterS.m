function y=filterS(t, x, symm)
    tlen = length(t); N0 = (tlen - 1)/2;
    N = size(x, 1);
    n = size(x, 2);
    
    if symm
        y = [x((N0+1):(-1):2, :) ; x; x((N-1):(-1):(N - N0), :)];
    else
        y = [x((N - N0 + 1):N, :); x; x(1:N0, :)];
    end
    for k=1:n 
        z = conv(t, y(:,k));
        x(:,k) = z((2*N0+1):(length(z)-2*N0));
    end
    y = x;