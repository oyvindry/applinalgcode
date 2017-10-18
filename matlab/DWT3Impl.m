function x=DWT3Impl(x, nres, wave_name, bd_mode, dual, transpose)

    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('transpose','var')) transpose = 0; end

    f = find_kernel(wave_name, 1, dual, transpose);

    if transpose
        x = IDWT3Impl_internal(x, nres, f, bd_mode);
    else
        x = DWT3Impl_internal(x, nres, f, bd_mode);
    end
end
