function x=dwt_impl(x, nres, wave_name, forward, dim, bd_mode, dual, transpose)
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % nres:      Number of resolutions.
    % wave_name: Name of the wavelet.
    % forward:   Whether to apply the forward or reverse transform. Default: 1
    % dim:       The dimension of the transform (1 for sound, 2 for images). Default: 1
    % bd_mode:   Boundary extension mode. Default: 1
    % dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: 0
    % transpose: Whether the transpose is to be taken. Default: 0
    
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('forward','var')) forward = 1; end
    if (~exist('dim','var')) dim  = 1; end
    if (~exist('transpose','var')) transpose = 0; end
    
    f=find_kernel(wave_name, forward, dual, transpose)
    if transpose
        forward = ~forward;
    end
    if forward
        if dim == 2
            x = DWT2Impl_internal(x, nres, f, bd_mode);
        elseif dim == 1
            x = DWTImpl_internal(x, nres, f, bd_mode);
        end 
    else
        if dim ==2
            x = IDWT2Impl_internal(x, nres, f, bd_mode);
        elseif dim == 1
            x = IDWTImpl_internal(x, nres, f, bd_mode);
        end
    end
end