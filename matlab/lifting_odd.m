function x=lifting_odd(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    x(2:2:(N-1), :) = lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if strcmpi(bd_mode, 'per')
        x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(N, :) = x(N, :) + lambda2*x(N-1, :);
    end
end

