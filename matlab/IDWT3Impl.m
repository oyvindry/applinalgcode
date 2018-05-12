function x=IDWT3Impl(x, m, wave_name, bd_mode, prefilter_mode, dual, transpose, data_layout)
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('prefilter_mode','var')) prefilter_mode = 'none'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    if (~exist('data_layout','var')) data_layout = 'resolution'; end
    
    [wav_propsx, dual_wav_propsx] = find_wav_props(m, wave_name, bd_mode, size(x,1));
    [wav_propsy, dual_wav_propsy] = find_wav_props(m, wave_name, bd_mode, size(x,2));
    [wav_propsz, dual_wav_propsz] = find_wav_props(m, wave_name, bd_mode, size(x,3));
    [wav_propsx, fx, prefilterx] = find_kernel(wav_propsx, dual_wav_propsx, 0, dual, transpose, prefilter_mode);
    [wav_propsy, fy, prefiltery] = find_kernel(wav_propsy, dual_wav_propsy, 0, dual, transpose, prefilter_mode);
    [wav_propsz, fz, prefilterz] = find_kernel(wav_propsz, dual_wav_propsz, 0, dual, transpose, prefilter_mode);
    if transpose % if transpose, then f will we a dwt_kernel, 
        x =  DWT3Impl_internal(x, m, fx, fy, fz, bd_mode, prefilterx, prefiltery, prefilterz, wav_propsx, wav_propsy, wav_propsz, data_layout);
    else
        x = IDWT3Impl_internal(x, m, fx, fy, fz, bd_mode, prefilterx, prefiltery, prefilterz, wav_propsx, wav_propsy, wav_propsz, data_layout);
    end
end
