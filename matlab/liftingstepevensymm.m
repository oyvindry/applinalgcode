function x=liftingstepevensymm(lambda, x, bd_mode)
    N = size(x, 1);
    if ~bd_mode
        assert(mod(N,2) == 0)
    end
    if bd_mode
        x(1, :) = x(1, :) + 2*lambda*x(2, :); % Symmetric extension
    else
        x(1, :) = lambda*(x(2, :) + x(N, :)) + x(1, :);
    end
    x(3:2:(N-1), :) = x(3:2:(N-1), :) + lambda*(x(2:2:(N-2), :) + x(4:2:N, :)); % This saves one multiplication
    if mod(N,2) == 1 % last is odd
        x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
    end