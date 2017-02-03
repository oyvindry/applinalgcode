function X = DWT2Impl(X, nres, wave_name, bd_mode, dual)
    if (~exist('bd_mode')) bd_mode = 1; end
    if (~exist('dual')) dual  = 0; end
    
    f = findDWTKernel(wave_name);
    M = size(X, 1); N = size(X, 2); sz = size(X);
    M0 = size(X, 1); N0 = size(X, 2);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = 0:(nres - 1)
        sz2(1) = M; Y2 = zeros(sz2);
        sz1(1) = N; Y1 = zeros(sz1);
        if length(sz1)==1
            Y1=zeros(sz1, 1); Y2=zeros(sz2, 1);
        end
        for n = 1:2^res:N0
            Y2(:, :) = X(1:2^res:M0, n, :);
            X(1:2^res:M0, n, :) = f(Y2, bd_mode, dual);
        end
        for m = 1:2^res:M0
            Y1(:, :) = X(m, 1:2^res:N0, :);
            X(m, 1:2^res:N0, :) = f(Y1, bd_mode, dual);
        end
        M = ceil(M/2); N = ceil(N/2);
    end
    
    X = reorganize_coefficients2(X, nres, 1);   
