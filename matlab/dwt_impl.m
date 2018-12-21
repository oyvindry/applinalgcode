function x=dwt_impl(x, wave_name, m, bd_mode, prefilter_mode, dims, dual, transpose, data_layout)
    % Main function for computing the DWT of a given signal. Can be used for
    % all signals up to dimension 3.  The dimension of the data may be one
    % higher than the dimension of the transform, in which case the last
    % dimension is used for parallel computation.
    %
    % Note that this function computes all quantities needed from scratch in
    % order to compute the DWT for the wavelet in question.  This can be
    % time-consuming, and can be avoided by using the functions find_wav_props,
    % find_kernel,  the internal DWT functions dwt1_impl_internal,
    % dwt2_impl_internal, dwt3_impl_internal, as well as Matlabs persistence
    % functions. An example with minimum set of parameters is as follows:
    % 
    % [wav_props, dual_wav_props] = find_wav_props(wave_name);
    % save('wav_props.mat', 'wav_props', 'dual_wav_props');
    % ...
    % load('wav_props.mat');
    % [f, prefilter] = find_kernel(wav_props, dual_wav_props, 1);
    % x = dwt1_impl_internal(x, f);
    %     
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'spline53' - Spline 5/3 wavelet
    %            'splinex.x' - Spline wavelet with given number of vanishing 
    %                          moments for each filter
    %            'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Daubechies orthnormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    % m:         Number of resolutions. Default: 1.
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'none'   - Take no extra action at the boundaries
    %            'bd'     - Boundary wavelets
    % prefilter_mode: Possible modes are:
    %            'none' (default)
    %            'filter'
    %            'bd_pre' - Boundary wavelets with preconditioning
    % dims:      The number of dimensions to apply the DWT to. Always applied
    %            to the first dimensions. Default: max(dim(x)-1,1). This means 
    %            that sound with many channels, and images with many colour 
    %            components default to a one- and two-dimensional DWT, 
    %            respectively
    % dual:      Whether to apply the dual wavelet rather than the wavelet 
    %            itself. Default: 0
    % transpose: Whether the transpose is to be taken. Default: 0
    % data_layout: How data should be assembled. Possible modes are:
    %            'resolution': Lowest resolution first (default)
    %            'time': Sort according to time

    if (~exist('m','var')) m = 1; end
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('prefilter_mode','var')) prefilter_mode = 'none'; end
    if (~exist('dims','var')  || dims == 0) 
        dims = 1;
        if length(size(x)) > 1
            dims = length(size(x)) - 1; 
        end
    end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    if (~exist('data_layout','var')) data_layout = 'resolution'; end

    [wav_propsx, dual_wav_propsx] = find_wav_props(wave_name, m, bd_mode, size(x,1));
    [fx, prefilterx] = find_kernel(wav_propsx, dual_wav_propsx, 1, dual, transpose, prefilter_mode);
    offsets = [wav_propsx.offset_L wav_propsx.offset_R];
    if dims == 1
        if transpose % if transpose, then f will we an idwt_kernel, 
            x = idwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout);
        else
            x =  dwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout);
        end
    else
        [wav_propsy, dual_wav_propsy] = find_wav_props(wave_name, m, bd_mode, size(x,2));
        [fy, prefiltery] = find_kernel(wav_propsy, dual_wav_propsy, 1, dual, transpose, prefilter_mode);
        offsets = [offsets; wav_propsy.offset_L wav_propsy.offset_R];
        if dims == 2
            if transpose % if transpose, then f will we an idwt_kernel, 
                x = idwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout);
            else
                x = dwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout);
            end
        else
            [wav_propsz, dual_wav_propsz] = find_wav_props(wave_name, m, bd_mode, size(x,3));
            [fz, prefilterz] = find_kernel(wav_propsz, dual_wav_propsz, 1, dual, transpose, prefilter_mode);
            offsets = [offsets; wav_propsz.offset_L wav_propsz.offset_R];
            if dims == 3 % if not give error message
                if transpose % if transpose, then f will we an idwt_kernel, 
                    x = idwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout);
                else
                    x =  dwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout);
                end
            end
        end
    end         
end
