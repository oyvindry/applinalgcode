function x = DWT3Impl_internal(x, nres, f, bd_mode)
    M = size(x, 1); N = size(x, 2); K = size(x, 3);
    M0 = size(x, 1); N0 = size(x, 2); K0 = size(x, 3);
    sz = size(x);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];

    for res = 0:(nres - 1)

        sz1(1) = N; 
        sz2(1) = M; 

        Y1 = zeros(sz1);         
        Y2 = zeros(sz2);
        Y3 = zeros([K,1,M0]);

        for n = 1:2^res:N0
            Y2(:, :) = x(1:2^res:M0, n, :);
            x(1:2^res:M0, n, :) = f(Y2, bd_mode);
        end

        for m = 1:2^res:M0
            Y1(:, :) = x(m, 1:2^res:N0, :);
            x(m, 1:2^res:N0, :) = f(Y1, bd_mode);
        end

        for n = 1:2^res:N0
            Y3(:,1,:) = permute(x(:,n, 1:2^res:K0), [3,2,1]);
            x(:,n, 1:2^res:K0) = permute(f(Y3, bd_mode), [3,2,1]);
        end

        M = ceil(M/2); 
        N = ceil(N/2);
        K = ceil(K/2);
    
    end
    x = reorganize_coeffs3_forward(x, nres);   
end

function Y=reorganize_coeffs3_forward(X, nres)
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
    
    Y(1:lc1, 1:lc2, 1:lc3) = X(inds1, inds2, inds3);
    for res = nres:(-1):1
        inds1 = (2^(res - 1) + 1):2^res:M;
        inds2 = (2^(res - 1) + 1):2^res:N;
        inds3 = (2^(res - 1) + 1):2^res:K;
        lw1 = length(inds1);
        lw2 = length(inds2);
        lw3 = length(inds3);
        
        % [1,0,0]
        Y((lc1 + 1):(lc1 + lw1), 1:lc2, 1:lc3) = X(inds1, 1:2^res:N, 1:2^res:K);
        
        % [1,1,0]
        Y((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), 1:lc3) = X(inds1, inds2, 1:2^res:K);
        
        % [0,1,0]
        Y(1:lc1, (lc2 + 1):(lc2 + lw2), 1:lc3) = X(1:2^res:M, inds2, 1:2^res:K);
        
        % [0,0,1]
        Y(1:lc1, 1:lc2, (lc3 + 1):(lc3 + lw3)) = X(1:2^res:M, 1:2^res:N, inds3);
        
        % [1,0,1]
        Y((lc1 + 1):(lc1 + lw1), 1:lc2, (lc3 + 1):(lc3 + lw3)) = X(inds1, 1:2^res:N, inds3);
        
        % [1,1,1]
        Y((lc1 + 1):(lc1 + lw1), (lc2 + 1):(lc2 + lw2), (lc3 + 1):(lc3 + lw3)) = X(inds1, inds2, inds3);
        
        % [0,1,1]
        Y(1:lc1, (lc2 + 1):(lc2 + lw2), (lc3 + 1):(lc3 + lw3)) = X(1:2^res:M, inds2, inds3);

        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
        lc3 = lc3 + lw3;
    end
end
