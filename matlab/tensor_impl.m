function x = tensor_impl(x, S1, S2)
    sz = 1:(length(size(x)));
    x(1:end, :) = S1(x(1:end, :));
    x = permute(x,[2 1 sz(3:end)]);
    x(1:end, :) = S2(x(1:end, :));
    x = permute(x,[2 1 sz(3:end)]);
end