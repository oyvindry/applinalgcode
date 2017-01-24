function x = DWTKernelpwl2(x, symm, dual)
    if dual
        x = x/sqrt(2);
        x = liftingstepevensymm(0.5, x, symm);
        x = liftingstepoddsymm(-0.25, x, symm);
    else
        x = x*sqrt(2);
        x = liftingstepoddsymm(-0.5, x, symm);
        x = liftingstepevensymm(0.25, x, symm);
    end