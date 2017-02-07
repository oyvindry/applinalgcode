function x=IDWT2Impl_internal(x, nres, f, bd_mode)
    x = reorganize_coeffs2_reverse(x, nres);   
    
    M = size(x, 1); N = size(x, 2); sz = size(x);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = (nres - 1):(-1):0
        sz1(1) = length(1:2^res:N); sz2(1) = length(1:2^res:M);
        Y1 = zeros(sz1); Y2 = zeros(sz2);
        if length(sz1)==1
            Y1=zeros(sz1, 1); Y2=zeros(sz2, 1);
        end
        for n = 1:2^res:N
            Y2(:, :) = x(1:2^res:M, n, :);
            x(1:2^res:M, n, :) = f(Y2(:, :), bd_mode);
        end
        for m = 1:2^res:M
            Y1(:, :) = x(m, 1:2^res:N, :);
            x(m, 1:2^res:N, :) = f(Y1(:, :), bd_mode);
        end
    end
end

function Y=reorganize_coeffs2_reverse(X, nres)
    M = size(X, 1);
    N = size(X, 2);
    Y = zeros(size(X));
    inds1 = 1:2^nres:M;
    inds2 = 1:2^nres:N;
    lc1 = length(inds1);
    lc2 = length(inds2);
    Y(inds1, inds2, :) = X(1:lc1, 1:lc2, :);
    for res = nres:(-1):1
        inds1 = (2^(res - 1) + 1):2^res:M;
        inds2 = (2^(res - 1) + 1):2^res:N;
        lw1 = length(inds1);
        lw2 = length(inds2);
        
        Y(inds1, 1:2^res:M, :) = X((lc1 + 1):(lc1 + lw1), 1:lc2, :);
        Y(inds1, (1+2^(res-1)):2^res:M, :) = X((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :);
        Y(1:2^res:M, inds2, :) = X(1:lc1, (lc2 + 1):(lc2 + lw2), :);
        
        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
    end    
end