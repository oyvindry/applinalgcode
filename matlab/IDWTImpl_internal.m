function x = IDWTImpl_internal(x, nres, f, bd_mode)
    x = reorganize_coeffs_reverse(x, nres);
    for res = (nres - 1):(-1):0
        x(1:2^res:end, :) = f(x(1:2^res:end, :), bd_mode);
    end
end

function y = reorganize_coeffs_reverse(x, nres)
    N = size(x,1);
    y = zeros(size(x));
    inds = 1:2^nres:N;
    lc = length(inds);
    y(inds, :) = x(1:lc, :);
    for res = nres:(-1):1
        inds = (2^(res - 1) + 1):2^res:N;
        lw = length(inds);
        y(inds, :) = x((lc + 1):(lc + lw), :);
        lc = lc + lw;
    end
end