function X=IDWT2Impl(X, nres, f, mode, dualarg)
    f = findIDWTKernel(wave_name);
    symm = 1;
    if nargin >= 4
        symm = mode;
    end
    dual = 0;
    if nargin >= 5
        dual = dualarg;
    end
    
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
            X(1:2^res:M, n, :) = f(Y2(:, :), symm, dual);
        end
        for m = 1:2^res:M
            Y1(:, :) = X(m, 1:2^res:N, :);
            X(m, 1:2^res:N, :) = f(Y1(:, :), symm, dual);
        end
    end
