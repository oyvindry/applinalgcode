function y=filter_impl(t, x, bd_mode)
    tlen = length(t); N0 = (tlen - 1)/2;
    N = size(x, 1);
    n = size(x, 2);
    
    if strcmpi(bd_mode, 'symm')
        y = [x((N0+1):(-1):2, :) ; x; x((N-1):(-1):(N - N0), :)];
    elseif strcmpi(bd_mode, 'per')
        y = [x((N - N0 + 1):N, :); x; x(1:N0, :)];
    elseif strcmpi(bd_mode, 'none') || strcmpi(bd_mode, 'bd')
        y = [ zeros(N0, n); x; zeros(N0, n)];
    end
    for k=1:n 
        z = conv(t, y(:,k));
        x(:,k) = z((2*N0+1):(length(z)-2*N0));
    end
    y = x;