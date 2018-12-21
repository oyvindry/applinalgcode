function x=tensor2_impl(x, fx, fy, bd_mode)
    sz = 1:(length(size(x)));
    x = fx(x, bd_mode);
    x = permute(x,[2 1 sz(3:end)]);
    x = fy(x, bd_mode);
    x = permute(x,[2 1 sz(3:end)]);
end