function x = DWTImpl_internal(x, nres, f, bd_mode)
    for res=0:(nres - 1)
        x(1:2^res:end, :) = f(x(1:2^res:end, :), bd_mode);
    end
    x = reorganize_coeffs_forward(x, nres);
end

function y = reorganize_coeffs_forward(x, nres)
    N = size(x,1);
    y = zeros(size(x));
    inds = 1:2^nres:N;
    lc = length(inds);
    y(1:lc, :) = x(inds, :);
    for res = nres:(-1):1
        inds = (2^(res - 1) + 1):2^res:N;
        lw = length(inds);
        y((lc + 1):(lc + lw), :) = x(inds, :);
        lc = lc + lw;
    end
end