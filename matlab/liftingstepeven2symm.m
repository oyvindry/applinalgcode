function x = liftingstepeven2symm(lambda1, lambda2, x, symm)
    N = size(x, 1);
    if ~symm
        assert(mod(N,2) == 0)
    end
    if symm
        x(1, :) =lambda1*2*x(2, :) + x(1, :) + lambda2*2*x(4, :); % Symmetric extension
        x(3, :) = lambda1*(x(4, :) + x(2, :)) + x(3, :) + lambda2*(x(6, :) + x(2, :));
        if mod(N,2)==1
            x(N, :)   = lambda1*2*x(N-1, :) + x(N, :) + lambda2*2*x(N-3, :);
            x(N-2, :) = lambda1*(x(N-3, :) + x(N-1, :)) + x(N-2, :) + lambda2*(x(N-5, :) + x(N-1, :));
        else
            x(N-1, :) = lambda1*(x(N-2, :) + x(N, :)) + x(N-1, :) + lambda2*(x(N-4, :) + x(N-2, :));
        end
        x(5:2:(N-3), :) = lambda1*(x(4:2:(N-4), :) + x(6:2:(N-2), :)) + x(5:2:(N-3), :) + lambda2*(x(2:2:(N-6), :) + x(8:2:N, :));
    else % N must be even
        x(1, :) = lambda1*(x(2, :) + x(N, :)) + x(1, :) + lambda2*(x(4, :) + x(N-2, :));
        x(3, :) = lambda1*(x(4, :) + x(2, :)) + x(3, :) + lambda2*(x(6, :) + x(N, :));
        x(N-1, :) = lambda1*(x(N-2, :) + x(N, :)) + x(N-1, :) + lambda2*(x(N-4, :) + x(2, :));
        x(5:2:(N-3), :) = lambda1*(x(4:2:(N-4), :) + x(6:2:(N-2), :)) + x(5:2:(N-3), :) + lambda2*(x(2:2:(N-6), :) + x(8:2:N, :));
    end