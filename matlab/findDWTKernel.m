function f= findDWTKernel(wave_name)
    % Find the DWTKernel corresponding to the given wavelet name 
    if strcmpi(wave_name, 'cdf97')
        f = @(x, mode, dual) DWTKernel97(x, mode, dual);
    elseif strcmpi(wave_name, 'cdf53')
        f = @(x, mode, dual) DWTKernel53(x, mode, dual);
    elseif strcmpi(wave_name, 'pwl0')
        f = @(x, mode, dual) DWTKernelpwl0(x, mode, dual);
    elseif strcmpi(wave_name, 'pwl2')
        f = @(x, mode, dual) DWTKernelpwl2(x, mode, dual);
    elseif strcmpi(wave_name, 'Haar')
        f = @(x, mode, dual) DWTKernelHaar(x, mode, dual);
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        filters = liftingfactortho(vm, 0, 0);
        f = @(x, mode, dual) DWTKernelOrtho(x, filters, mode, dual);
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, mode, dual) DWTKernelOrtho(x, filters, mode, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        filters = liftingfactortho(vm, 1, 0);
        f = @(x, mode, dual) DWTKernelOrtho(x, filters, mode, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, mode, dual) DWTKernelOrtho(x, filters, mode, dual);
    end

