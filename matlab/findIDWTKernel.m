function f = findIDWTKernel(wave_name)
    % Find the IDWTKernel corresponding to the given wavelet name 
    if strcmpi(wave_name, 'cdf97')
        f = @IDWTKernel97;
    elseif strcmpi(wave_name, 'cdf53')
        f = @IDWTKernel53;
    elseif strcmpi(wave_name, 'pwl0')
        f = @IDWTKernelpwl0;
    elseif strcmpi(wave_name, 'pwl2')
        f = @IDWTKernelpwl2;
    elseif strcmpi(wave_name, 'Haar')
        f = @IDWTKernelHaar;
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        filters = getDBfilter(vm, 0);
        f = @(x, bd_mode, dual) IDWTKernelOrtho(x, filters, bd_mode, dual);
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, bd_mode, dual) IDWTKernelOrtho(x, filters, bd_mode, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        filters = getDBfilter(vm, 1);
        f = @(x, bd_mode, dual) IDWTKernelOrtho(x, filters, bd_mode, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, bd_mode, dual) IDWTKernelOrtho(x, filters, bd_mode, dual);
    end
end
