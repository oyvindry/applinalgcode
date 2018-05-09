function x=DWTImpl(x, m, wave_name,  bd_mode, prefilter_mode, dual, transpose, data_layout)
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % m:         Number of resolutions.
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'cdf53' - Spline 5/3 wavelet
    %            'splinex.x' - Spline wavelet with given number of vanishing moments for each filter
    %            'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Daubechies orthnormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'none'   - Take no extra action at the boundaries
    %            'bd'     - Boundary wavelets
    % prefilter_mode: Possible modes are:
    %            'none' (default)
    %            'filter'
    %            'bd_pre' - Boundary wavelets with preconditioning
    % dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: 0
    % transpose: Whether the transpose is to be taken. Default: 0
    % data_layout: How data should be assembled. Possible modes are:
    %            'resolution': Lowest resolution first (default)
    %            'time': Sort according to time

    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('prefilter_mode','var')) prefilter_mode = 'none'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    if (~exist('data_layout','var')) data_layout = 'resolution'; end
    
    [wav_props, dual_wav_props] = find_wav_props(m, wave_name, bd_mode, size(x,1));
    [wav_props, f, prefilter] = find_kernel(wav_props, dual_wav_props, 1, dual, transpose, prefilter_mode);
    if transpose % if transpose, then f will we an idwt_kernel, 
        x = IDWTImpl_internal(x, m, f, bd_mode, prefilter, wav_props, data_layout);     
    else
        x = DWTImpl_internal(x, m, f, bd_mode, prefilter, wav_props, data_layout);
    end
end
