function x=tensor3_impl(x, fx, fy, fz, bd_mode)
    sz = 1:(length(size(x)));
    x = fx(x, bd_mode);
    x = permute(x,[2 3 1 sz(4:end)]);
    x = fy(x, bd_mode);
    x = permute(x,[2 3 1 sz(4:end)]);
    x = fz(x, bd_mode);
    x = permute(x,[2 3 1 sz(4:end)]);
end