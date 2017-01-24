function y = reorganize_coefficients(x, nres, forward)
    N = size(x,1);
    y = zeros(size(x));
    inds = 1:2^nres:N;
    lc = length(inds);
    if forward
        y(1:lc, :) = x(inds, :);
    else
        y(inds, :) = x(1:lc, :);
    end
    for res = nres:(-1):1
        inds = (2^(res - 1) + 1):2^res:N;
        lw = length(inds);
        if forward
            y((lc + 1):(lc + lw), :) = x(inds, :);
        else
            y(inds, :) = x((lc + 1):(lc + lw), :);
        end
        lc = lc + lw;
    end