function x=DWTKernelFilters(H0, H1, G0, G1, x, symm, dual)
    f0 = H0; f1 = H1;
    if dual
        f0 = G0; f1 = G1;
    end  
    N = length(x);
    x0 = filterS(f0, x, symm);
    x1 = filterS(f1, x, symm);
    x(1:2:N) = x0(1:2:N);
    x(2:2:N) = x1(2:2:N);
