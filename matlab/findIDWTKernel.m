function f = findIDWTKernel(wave_name)
    % Find the IDWTKernel corresponding to the given wavelet name 
    if strcmpi(wave_name, 'cdf97')
        f = @(x, symm, dual) IDWTKernel97(x, symm, dual);
    elseif strcmpi(wave_name, 'cdf53')
        f = @(x, symm, dual) IDWTKernel53(x, symm, dual);
    elseif strcmpi(wave_name, 'pwl0')
        f = @(x, symm, dual) IDWTKernelpwl0(x, symm, dual);
    elseif strcmpi(wave_name, 'pwl2')
        f = @(x, symm, dual) IDWTKernelpwl2(x, symm, dual);
    elseif strcmpi(wave_name, 'Haar')
        f = @(x, symm, dual) IDWTKernelHaar(x, symm, dual);
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        filters = liftingfactortho(vm, 0, 0);
        f = @(x, symm, dual) IDWTKernelOrtho(x, filters, symm, dual);
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, symm, dual) IDWTKernelOrtho(x, filters, symm, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        filters = liftingfactortho(vm, 1, 0);
        f = @(x, symm, dual) IDWTKernelOrtho(x, filters, symm, dual);
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, symm, dual) IDWTKernelOrtho(x, filters, symm, dual);
    end
end
