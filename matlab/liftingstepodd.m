function x = liftingstepodd(lambda1, lambda2, x, mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    x(2:2:(N-1), :) = lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if mode==0
        x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
    elseif mode>=2
        x(N, :) = x(N, :) + lambda2*x(N-1, :);
    end
  
