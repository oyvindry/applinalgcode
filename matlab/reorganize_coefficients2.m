function Y=reorganize_coefficients2(X, nres, forward)
    M = size(X, 1);
    N = size(X, 2);
    Y = zeros(size(X));
    inds1 = 1:2^nres:M;
    inds2 = 1:2^nres:N;
    lc1 = length(inds1);
    lc2 = length(inds2);
    if forward
        Y(1:lc1, 1:lc2, :) = X(inds1, inds2, :);
    else
        Y(inds1, inds2, :) = X(1:lc1, 1:lc2, :);
    end
    for res = nres:(-1):1
        inds1 = (2^(res - 1) + 1):2^res:M;
        inds2 = (2^(res - 1) + 1):2^res:N;
        lw1 = length(inds1);
        lw2 = length(inds2);
        if forward
            Y((lc1 + 1):(lc1 + lw1), 1:lc2, :) = X(inds1, 1:2^res:M, :);
            Y((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :) = X(inds1, (1+2^(res-1)):2^res:M, :);
            Y(1:lc1, (lc2 + 1):(lc2 + lw2), :) = X(1:2^res:M, inds2, :);
        else
            Y(inds1, 1:2^res:M, :) = X((lc1 + 1):(lc1 + lw1), 1:lc2, :);
            Y(inds1, (1+2^(res-1)):2^res:M, :) = X((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :);
            Y(1:2^res:M, inds2, :) = X(1:lc1, (lc2 + 1):(lc2 + lw2), :);
        end
        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
    end