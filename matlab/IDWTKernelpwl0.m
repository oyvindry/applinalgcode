function x = IDWTKernelpwl0(x, symm, dual)
    if dual
        x = x*sqrt(2);
        x = liftingstepevensymm(-0.5, x, symm);
    else
        x = x/sqrt(2);
        x = liftingstepoddsymm(0.5, x, symm);
    end