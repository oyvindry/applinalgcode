function x = DWTKernel97(x,symm,dual)
    lambda1 = -0.586134342059950;
    lambda2 = -0.668067171029734;
    lambda3 = 0.070018009414994;
    lambda4 = 1.200171016244178;
    alpha = -1.149604398860250;
    beta = -0.869864451624777;
    N = size(x, 1);
  
    if dual
        x(1:2:N, :) =x(1:2:N, :)/alpha;
        x(2:2:N, :)=x(2:2:N, :)/beta;
        x = liftingstepevensymm(lambda4, x, symm);
        x = liftingstepoddsymm(lambda3, x, symm);
        x = liftingstepevensymm(lambda2, x, symm);
        x = liftingstepoddsymm(lambda1, x, symm);
    else
        x(1:2:N, :) = x(1:2:N, :)*alpha;
        x(2:2:N, :) = x(2:2:N, :)*beta;
        x = liftingstepoddsymm(-lambda4, x, symm);
        x = liftingstepevensymm(-lambda3, x, symm);
        x = liftingstepoddsymm(-lambda2, x, symm);
        x = liftingstepevensymm(-lambda1, x, symm);
    end