function x = liftingstepodd2symm(lambda1, lambda2, x, symm)
    N = size(x, 1);
    if ~symm
        assert(mod(N,2) == 0)
    end
    if symm
        x(2, :) = lambda1*2*(x(1, :) + x(3, :)) + x(2, :) + lambda2*(x(5, :) + x(3, :)); % Symmetric extension
        if mod(N,2) == 1
            x(N-1, :) = lambda1*(x(N-2, :) + x(N, :)) + x(N-1, :) + lambda2*(x(N-4, :) + x(N-2, :)); 
        else
        x(N, :)   = lambda1*2*x(N-1, :) + x(N, :) + lambda2*2*x(N-3, :);
        x(N-2, :) = lambda1*(x(N-3, :) + x(N-1, :)) + x(N-2, :) + lambda2*(x(N-5, :) + x(N-1, :));
    end
    x(4:2:(N-3), :) = lambda1*(x(3:2:(N-4), :) + x(5:2:(N-2), :)) + x(4:2:(N-3), :) + lambda2*(x(1:2:(N-6), :) + x(7:2:N, :));
    else % N must be even
        x(2, :) = lambda1*(x(1, :) + x(3, :)) + x(2, :) + lambda2*(x(5, :) + x(3, :));
        x(N-2, :) = lambda1*(x(N-3, :) + x(N-1, :)) + x(N-2, :) + lambda2*(x(N-5, :) + x(N-1, :));
        x(N, :) = lambda1*2*x(N-1, :) + x(N, :) + lambda2*2*x(N-3, :);
        x(4:2:(N-4), :) = lambda1*(x(3:2:(N-5), :) + x(5:2:(N-3), :)) + x(4:2:(N-4), :) + lambda2*(x(1:2:(N-7), :) + x(7:2:(N-1), :));
    end