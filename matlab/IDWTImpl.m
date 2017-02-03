function x = IDWTImpl(x, nres, wave_name, mode, dual)
    % See DWTImpl
    if (~exist('mode')) mode = 1; end
    if (~exist('dual')) dual  = 0; end
    
    f = findIDWTKernel(wave_name);
    x = reorganize_coefficients(x, nres, 0);
    N = size(x, 1);
    for res = (nres - 1):(-1):0
        x(1:2^res:N, :) = f(x(1:2^res:N, :), mode, dual);
    end
end

