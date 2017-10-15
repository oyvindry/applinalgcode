function x=lifting_even(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    if strcmpi(bd_mode, 'per')
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
end

