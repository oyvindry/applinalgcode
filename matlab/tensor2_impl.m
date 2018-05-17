function x=tensor2_impl(x, indsx, indsy, fx, fy, bd_mode)
    sz = 1:(length(size(x)));
    x(indsx, :) = fx(x(indsx, :), bd_mode);
    x = permute(x,[2 1 sz(3:end)]);
    x(indsy, :) = fy(x(indsy, :), bd_mode);
    x = permute(x,[2 1 sz(3:end)]);
end