function x=IDWTKernelFilters(H0, H1, G0, G1, x, symm, dual)
    f0 = G0; f1 = G1;
    if dual
        f0 = H0; f1 = H1;
    end
    N = length(x);
    x0 = x; x0(2:2:N) = 0;
    x1 = x; x1(1:2:N) = 0;
    x0 = filterS(f0, x0, symm);
    x1 = filterS(f1, x1, symm);
    x = x0 + x1;