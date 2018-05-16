function X = tensor_impl(X, S1, S2)
    sz = 1:(length(size(X)));
    x(1:end, :) = S1(x(1:end, :), bd_mode);
    x = permute(x,[2 1 sz(3:end)]);
    x(1:end, :) = S2(x(1:end, :), bd_mode);
    x = permute(x,[2 1 sz(3:end)]);