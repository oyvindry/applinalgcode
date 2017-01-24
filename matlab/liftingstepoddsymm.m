function x=liftingstepoddsymm(lambda, x, symm)
    N = size(x, 1);
    if ~symm
        assert(mod(N,2) == 0)
    end
    x(2:2:(N-1), :) = x(2:2:(N-1), :) + lambda*(x(1:2:(N-2), :) + x(3:2:N, :)); % This saves one multiplication
    if mod(N,2)==0 % last is even
        if symm
            x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
        else
            x(N, :) = lambda*(x(1, :) + x(N-1, :)) + x(N, :);
        end
    end