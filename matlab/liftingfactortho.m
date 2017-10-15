function filters=liftingfactortho(N, type, debug_mode)

    % Computes the filter coefficients of orthonormal wavelets with N vanishing
    % moments.
    %
    % N:    Number of vanishing moments
    % type: The type of orthonormal wavelet.
    %       0: Daubechies wavelets with minimum phase (default)  
    %       1: Symlets - wavelets with close to linear phase (almost symetric) 
    % debug_mode: Wheter or not this function should be ran in debug mode
    % 
    if (nargin <  3)
        debug_mode = 0;
    end
    if (nargin == 1)
        type = 0;
    end


        
    % First the right edge
    % h0 - low pass filter
    % h1 - high pass filter
    if (type == 0)
        [h0, h1, g0, g1] = h0h1computeortho(N);
    elseif (type == 1)
        [h0, h1, g0, g1] = h0h1computesym(N);
    end
     
    h0 = flip(h0);
    h1 = flip(h1);
    g0 = flip(g0);
    g1 = flip(g1);
    %filters = liftingstepscomputeortho(h0, h1);
     
    [W, A_pre, A_pre_inv] = bw_compute_left(h0, g0, debug_mode); % Lower right (3N-1)x(2N) matrix
    
    %%%%%%%%%%%%%%%%%%%%%%% Save to file %%%%%%%%%%%%%%%%%%%%%%%%
    if (type == 0)
        dest = 'my_CDJV_filters';
        
        filenameW = sprintf('%s/WR_db%d.mat', dest, N);
        filenameA_pre = sprintf('%s/AR_pre_db%d.mat', dest, N);
        filenameA_pre_inv = sprintf('%s/AR_pre_inv_db%d.mat', dest, N);
        
        WR = W;
        AR_pre = A_pre;
        AR_pre_inv = A_pre_inv;
        save(filenameW, 'WR');
        save(filenameA_pre, 'AR_pre');
        save(filenameA_pre_inv, 'AR_pre_inv');
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WR = zeros(size(W));
    for k=1:N
        WR(:,[2*k-1 2*k]) = W(size(W,1):(-1):1,2*N+1-[2*k 2*k-1]); 
    end
     
    % Then the left edge
    h0 = flip(h0);
    h1 = flip(h1);
    g0 = flip(g0);
    g1 = flip(g1);
     
    filters = liftingstepscomputeortho(h0, h1);
    filters.A_R_pre = fliplr(flipud(A_pre));
    filters.A_R_pre_inv = fliplr(flipud(A_pre_inv));
    [WL,A_pre, A_pre_inv] = bw_compute_left(h0, g0, debug_mode); % Upper left (3N-1)x(2N) matrix
    if (type == 0)
        dest = 'my_CDJV_filters';
        
        filenameW = sprintf('%s/WL_db%d.mat', dest, N);
        filenameA_pre = sprintf('%s/AL_pre_db%d.mat', dest, N);
        filenameA_pre_inv = sprintf('%s/AL_pre_inv_db%d.mat', dest, N);
        
        AL_pre = A_pre;
        AL_pre_inv = A_pre_inv;
        save(filenameW, 'WL');
        save(filenameA_pre, 'AL_pre');
        save(filenameA_pre_inv, 'AL_pre_inv');
    end  
    
    filters.A_L_pre = A_pre;
    filters.A_L_pre_inv = A_pre_inv;
     
    % Compute the left and right parts of the IDWT for boundary handling
    M = 6*N;
    seg1 = zeros(M); % One bigger than is actually needed
     
    filters.AL = zeros(size(WL));
    filters.AR = zeros(size(WR));
    
    for k=0:(M-1)
        x = zeros(M,1);
        x(k+1) = 1;
        seg1(:,k+1) = idwt_kernel_ortho(x, filters, 'bd');
    end

    % Why is WL computed in another way than WR? 
    [w1, w2] = size(WL);
    filters.AL=WL-seg1(1:w1,1:w2);
    filters.AR=WR-seg1((M-w1+1):M,(M-w2+1):M);

end


