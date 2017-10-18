function x=IDW2TImpl(x, nres, wave_name, bd_mode, dual, transpose)
    % x:         Matrix whose IDWT will be computed along the first dimension(s).      
    % nres:      Number of resolutions.
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'cdf53' - Spline 5/3 wavelet  
    %            'splinex.x' - Spline wavelet with given number of vanishing moments for each filter
    %            'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Dauberchies orthnormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'bd'     - Boundary wavelets
    %            'bd_pre' - Boundary wavelets with preconditioning
    % dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: 0
    % transpose: Whether the transpose is to be taken. Default: 0
    
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    
    f = find_kernel(wave_name, 0, dual, transpose);
    
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
        if (strcmpi(bd_mode, 'bd_pre')) % Apply preconditioning
            n = size(filters.A_L_pre,1);
            x(1:n, :) = filters.A_L_pre*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre*x((end-n+1):end, :);
            
            x(:, 1:n) = (filters.A_L_pre*x(:, 1:n)')';
            x(:, (end-n+1):end) = (filters.A_R_pre*x(:, (end-n+1):end)')';
        end
        x = DWT2Impl_internal(x, nres, f, bd_mode);
    else
        x = IDWT2Impl_internal(x, nres, f, bd_mode);
        
        if (strcmpi(bd_mode, 'bd_pre')) % Apply postconditioning
            n = size(filters.A_L_pre_inv,1);
            
            x(1:n, :) = filters.A_L_pre_inv*x(1:n, :);
            x((end-n+1):end, :) = filters.A_R_pre_inv*x((end-n+1):end, :);
            
            x(:, 1:n) = (filters.A_L_pre_inv*(x(:, 1:n)'))';
            x(:, (end-n+1):end) = (filters.A_R_pre_inv*(x(:, (end-n+1):end)'))';
        end
    end
end
