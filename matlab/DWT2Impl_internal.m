function x = DWT2Impl_internal(x, nres, f, bd_mode)
    M = size(x, 1); N = size(x, 2); sz = size(x);
    M0 = size(x, 1); N0 = size(x, 2);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = 0:(nres - 1)
        
        sz2(1) = M; 
        sz1(1) = N; 
        
        if length(sz1)==1
            Y1=zeros(sz1, 1); 
            Y2=zeros(sz2, 1);
        else
            Y2 = zeros(sz2);
            Y1 = zeros(sz1);         
        end
        
        for n = 1:2^res:N0
            Y2(:, :) = x(1:2^res:M0, n, :);
            x(1:2^res:M0, n, :) = f(Y2, bd_mode);
        end
        
        for m = 1:2^res:M0
            Y1(:, :) = x(m, 1:2^res:N0, :);
            x(m, 1:2^res:N0, :) = f(Y1, bd_mode);
        end
        
        M = ceil(M/2); N = ceil(N/2);
    end
    
    x = reorganize_coeffs2_forward(x, nres);   
end

function Y=reorganize_coeffs2_forward(X, nres)
    M = size(X, 1);
    N = size(X, 2);
    Y = zeros(size(X));
    inds1 = 1:2^nres:M;
    inds2 = 1:2^nres:N;
    lc1 = length(inds1);
    lc2 = length(inds2);
    
    Y(1:lc1, 1:lc2, :) = X(inds1, inds2, :);
    for res = nres:(-1):1
        inds1 = (2^(res - 1) + 1):2^res:M;
        inds2 = (2^(res - 1) + 1):2^res:N;
        lw1 = length(inds1);
        lw2 = length(inds2);
        
        Y((lc1 + 1):(lc1 + lw1), 1:lc2, :) = X(inds1, 1:2^res:N, :);
        Y((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :) = X(inds1, inds2, :);
        Y(1:lc1, (lc2 + 1):(lc2 + lw2), :) = X(1:2^res:M, inds2, :);

        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
    end    
end
