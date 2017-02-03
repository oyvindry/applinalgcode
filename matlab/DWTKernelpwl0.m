function x = DWTKernelpwl0(x, bd_mode, dual)
    if dual
        x = x/sqrt(2);
        x = liftingstepevensymm(0.5, x, bd_mode);
    else
        x = x*sqrt(2);
        x = liftingstepoddsymm(-0.5, x, bd_mode);
    end