function x=IDWT3Impl_internal(x, nres, f, bd_mode)
    x = reorganize_coeffs3_reverse(x, nres);   
    
    M = size(x, 1); N = size(x, 2); K = size(x,3);
    sz = size(x);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    sz3 = [1,1,M];
    for res = (nres - 1):(-1):0
        sz1(1) = length(1:2^res:N); 
        sz2(1) = length(1:2^res:M);
        sz3(1) = length(1:2^res:K);
        
        Y1 = zeros(sz1);
        Y2 = zeros(sz2);
        Y3 = zeros(sz3);

        for n = 1:2^res:N
            Y2(:, :) = x(1:2^res:M, n, :);
            x(1:2^res:M, n, :) = f(Y2(:, :), bd_mode);
        end

        for m = 1:2^res:M
            Y1(:, :) = x(m, 1:2^res:N, :);
            x(m, 1:2^res:N, :) = f(Y1(:, :), bd_mode);
        end
        
        for n = 1:2^res:N
            Y3(:,:) = permute(x(:, n, 1:2^res:K),[3,2,1]);
            x(:, n, 1:2^res:K) = permute(f(Y3(:,:), bd_mode), [3,2,1]);
        end
    end
end


function Y=reorganize_coeffs3_reverse(X, nres)
    M = size(X, 1);
    N = size(X, 2);
    K = size(X, 3);
    Y = zeros(size(X));
    inds1 = 1:2^nres:M;
    inds2 = 1:2^nres:N;
    inds3 = 1:2^nres:K;
    lc1 = length(inds1);
    lc2 = length(inds2);
    lc3 = length(inds3);
    Y(inds1, inds2, inds3) = X(1:lc1, 1:lc2, 1:lc3);
    for res = nres:(-1):1

        inds1 = (2^(res-1) + 1):2^res:M;
        inds2 = (2^(res-1) + 1):2^res:N;
        inds3 = (2^(res-1) + 1):2^res:K;
        lw1 = length(inds1);
        lw2 = length(inds2);
        lw3 = length(inds3);

        Y(inds1, 1:2^res:N, inds3) = X((lc1 + 1):(lc1 + lw1), 1:lc2, (lc3 + 1):(lc3 + lw3));
        Y(inds1, inds2, inds3)     = X((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), (lc3 + 1):(lc3 + lw3));
        Y(1:2^res:M, inds2, inds3) = X(1:lc1, (lc2 + 1):(lc2 + lw2), (lc3 + 1):(lc3 + lw3));
        
        Y(1:2^res:M, 1:2^res:K, inds3) = X(1:lc1, 1:lc2, (lc3 + 1):(lc3 + lw3));

        Y(inds1, 1:2^res:N, 1:2^res:K) = X((lc1 + 1):(lc1 + lw1), 1:lc2, 1:lc3);
        Y(inds1, inds2, 1:2^res:K)     = X((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), 1:lc3);
        Y(1:2^res:M, inds2, 1:2^res:K) = X(1:lc1, (lc2 + 1):(lc2 + lw2), 1:lc3);

        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
        lc3 = lc3 + lw3;

    end
end




