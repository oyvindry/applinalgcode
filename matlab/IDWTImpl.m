function x = IDWTImpl(x, nres, wave_name, mode, dualarg)
    % See DWTImpl
    if (~exist('mode')) mode = 1; end
    if (~exist('dualarg')) dualarg  = 0; end
    
    f = findIDWTKernel(wave_name);
    x = IDWTReadKernel(x, nres, f, mode, dualarg);
end

