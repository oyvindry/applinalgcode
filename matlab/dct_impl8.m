function x = dct_impl8(x, bd_mode)
    N = size(x, 1);
    for n = 1:8:N
        x(n:(n+7), :) = dct(x(n:(n+7), :));
    end
end