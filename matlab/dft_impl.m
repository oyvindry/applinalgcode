function y = dft_impl(x, forward)
    sign = -1;
    if nargin >= 2 & forward == 0 
        sign = 1;
    end
    N = size(x, 1);
    y = zeros(size(x));
    for n = 1:N
        D = exp(-sign*2*pi*1i*(n-1)*(0:(N-1))/N);
        for s2 = 1:size(x,2)
            y(n, s2) = dot(D, x(:, s2));
        end
    end
    if sign == 1
        y = y/N;
    end
end