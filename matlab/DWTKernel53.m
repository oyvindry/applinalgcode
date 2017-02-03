function x = DWTKernell53(x, bd_mode, dual)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
    N = size(x, 1);
  
    if dual
        x(1:2:N, :) = x(1:2:N, :)/alpha;
        x(2:2:N, :) = x(2:2:N, :)/beta;
        x = liftingstepevensymm(lambda2, x, bd_mode);
        x = liftingstepoddsymm(lambda1, x, bd_mode);
    else
        x(1:2:N, :) = x(1:2:N, :)*alpha;
        x(2:2:N, :) = x(2:2:N, :)*beta;
        x = liftingstepoddsymm(-lambda2, x, bd_mode);
        x = liftingstepevensymm(-lambda1, x, bd_mode);
    end