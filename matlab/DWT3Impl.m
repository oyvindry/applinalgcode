function x=DWT3Impl(x, nres, wave_name, bd_mode, dual, transpose)

    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end

    f = find_kernel(wave_name, 1, dual, transpose);
    
    if (strcmpi(bd_mode, 'bd_pre')) % Apply postconditioning
        if (strcmpi(wave_name(1:2),'db') )
            vm = str2num(wave_name(3:end));
            filters = getDBfilter(vm, 0);
        end
        
        if( strcmpi(wave_name(1:3),'sym') )
            vm = str2num(wave_name(4:end));
            filters = getDBfilter(vm, 0);
        end
    end

    if transpose
        x = IDWT3Impl_internal(x, nres, f, bd_mode);
        
        if (strcmpi(bd_mode, 'bd_pre')) % Apply postconditioning
            n = size(filters.A_L_pre_inv,1);
            x(1:n, :) = filters.A_L_pre_inv*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre_inv*x((end-n+1):end, :);
            
            x = permute(x, [2,1,3]);
            x(1:n, :) = filters.A_L_pre_inv*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre_inv*x((end-n+1):end, :);
            
            x = permute(x, [3,2,1]);
            x(1:n, :) = filters.A_L_pre_inv*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre_inv*x((end-n+1):end, :);
            
            x = permute(x, [2,3,1]);
        end
    else
        if (strcmpi(bd_mode, 'bd_pre')) % Apply preconditioning
            n = size(filters.A_L_pre,1);
            x(1:n, :) = filters.A_L_pre*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre*x((end-n+1):end, :);
            x = permute(x, [2,1,3]);
            x(1:n, :) = filters.A_L_pre*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre*x((end-n+1):end, :);
            x = permute(x, [3,2,1]);
            x(1:n, :) = filters.A_L_pre*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre*x((end-n+1):end, :);
            x = permute(x, [2,3,1]);
        end
        
        x = DWT3Impl_internal(x, nres, f, bd_mode);
    end
end
