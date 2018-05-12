function x=DWT2Impl(x, m, wave_name, bd_mode, prefilter_mode, dual, transpose, data_layout)
    
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('prefilter_mode','var')) prefilter_mode = 'none'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    if (~exist('data_layout','var')) data_layout = 'resolution'; end


    [wav_propsx, dual_wav_propsx] = find_wav_props(m, wave_name, bd_mode, size(x,1));
    [wav_propsy, dual_wav_propsy] = find_wav_props(m, wave_name, bd_mode, size(x,2));
    [wav_propsx, fx, prefilterx] = find_kernel(wav_propsx, dual_wav_propsx, 1, dual, transpose, prefilter_mode);
    [wav_propsy, fy, prefiltery] = find_kernel(wav_propsy, dual_wav_propsy, 1, dual, transpose, prefilter_mode);
    if transpose % if transpose, then f will we an idwt_kernel, 
        x = IDWT2Impl_internal(x, m, fx, fy, bd_mode, prefilterx, prefiltery, wav_propsx, wavpropsy, data_layout);     
    else
        x =  DWT2Impl_internal(x, m, fx, fy, bd_mode, prefilterx, prefiltery, wav_propsx, wavpropsy, data_layout);
    end
end