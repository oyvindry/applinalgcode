function X=IDWT2Impl(X, nres, wave_name, bd_mode, dual)
    if (~exist('bd_mode')) bd_mode = 1; end
    if (~exist('dual')) dual  = 0; end
    
    f = findIDWTKernel(wave_name);
    X = reorganize_coefficients2(X, nres, 0);   
    
    M = size(X, 1); N = size(X, 2); sz = size(X);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = (nres - 1):(-1):0
        sz1(1) = length(1:2^res:N); sz2(1) = length(1:2^res:M);
        Y1 = zeros(sz1); Y2 = zeros(sz2);
        if length(sz1)==1
            Y1=zeros(sz1, 1); Y2=zeros(sz2, 1);
        end
        for n = 1:2^res:N
            Y2(:, :) = X(1:2^res:M, n, :);
            X(1:2^res:M, n, :) = f(Y2(:, :), bd_mode, dual);
        end
        for m = 1:2^res:M
            Y1(:, :) = X(m, 1:2^res:N, :);
            X(m, 1:2^res:N, :) = f(Y1(:, :), bd_mode, dual);
        end
    end
