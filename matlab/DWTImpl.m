function x=DWTImpl(x, nres, wave_name, bd_mode, dual, transpose)
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % nres:      Number of resolutions.
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'cdf53' - Spline 5/3 wavelet  
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
    
    f = find_kernel(wave_name, 1, dual, transpose);
    if transpose
        x = IDWTImpl_internal(x, nres, f, bd_mode);
    else
        x = DWTImpl_internal(x, nres, f, bd_mode);
    end
end
