function X = tensor_impl(X, S1, S2)
    M = size(X, 1); N = size(X, 2); sz = size(X);
    sz1 = sz; sz1(1) = []; sz1 = [sz1 1];
    sz2 = sz; sz2(2) = []; sz2 = [sz2 1];
    Y1 = zeros(sz1);
    Y2 = zeros(sz2);
    
    for n = 1:N
        Y2(:, :) = X(:, n, :);
        X(:, n, :) = S1(Y2);
    end
    for m = 1:M
        Y1(:, :) = X(m, :, :);
        X(m, :, :) = S2(Y1);
    end