function x = liftingstepeven(lambda1, lambda2, x, mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    if mode==0
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif mode>=2
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
  
