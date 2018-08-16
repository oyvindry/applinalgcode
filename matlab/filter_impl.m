function x=filter_impl(t, x, bd_mode)
    N0 = (length(t) - 1)/2;
    N = size(x, 1);
    szx = size(x);
    n = prod(szx(2:end));
    
    if strcmpi(bd_mode, 'symm')
        y = [x((N0+1):(-1):2, :) ; x(:,:); x((N-1):(-1):(N - N0), :)];
    elseif strcmpi(bd_mode, 'per')
        y = [x((N - N0 + 1):N, :); x(:,:); x(1:N0, :)];
    elseif strcmpi(bd_mode, 'none') || strcmpi(bd_mode, 'bd')
        y = [ zeros(N0, n); x(:, :); zeros(N0, n)];
    end
    for k=1:size(y,2)
        z = conv(t, y(:,k));
        x(:,k) = z((2*N0+1):(length(z)-2*N0));
    end
end