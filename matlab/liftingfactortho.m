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
    
    % We remove the persistent variables until we are done testing, so that
    % everything is recomputed each time one call this function.

    %persistent filterMap;
    %if (isempty(filterMap)) 
    %    filterMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    %end
    %if (filterMap.isKey(N) && debugMode == 0) 
    %    filters = filterMap(N);
    %else
        
    % First the right edge
    if (type == 0)
        [h0, h1, g0, g1] = h0h1computeortho(N);
    elseif (type == 1)
        [h0, h1, g0, g1] = h0h1computesym(N);
    end
    h0 = flip(h0);
    h1 = flip(h1);
    g0 = flip(g0);
    g1 = flip(g1);
    filters = liftingstepscomputeortho(h0, h1);
    
    [W, A_pre, A_pre_inv] = bw_compute_left(h0, g0, debug_mode); % Lower right (3N-1)x(2N) matrix
    %filters.A_R_pre = fliplr(flipud(A_pre));
    %filters.A_R_pre_inv = fliplr(flipud(A_pre_inv));
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
        seg1(:,k+1) = IDWTKernelOrtho(x, filters, 2, 0);
    end
    
    [w1, w2] = size(WL);
    filters.AL=WL-seg1(1:w1,1:w2);
    filters.AR=WR-seg1((M-w1+1):M,(M-w2+1):M);
    
        %if (debugMode == 0)
        %    % Store filters to current session
        %    filterMap(N) = filters;
        %end
    %end
end

function [h0, h1, g0, g1]=h0h1computeortho(N)
    % Comptues the wavelet coefficients of the orthonormal Daubechies wavelet
    % N vanishing moments and with minimum phase   
    vals=computeQN(N);
    rts=roots(vals)';
    rts1=rts(find(abs(rts)>1));

    g0=1;
    for rt=rts1
        g0=conv(g0,[-rt 1]);
    end
    g0 = real(g0);
    K=sqrt(vals(1)*(-1)^(length(rts1))/abs(prod(rts1)));
    g0=K*g0;
    for k=1:N
        g0=conv(g0,[1/2 1/2]);
    end
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); 
    h1=fliplr(g1);
end

function [h0, h1, g0, g1]=h0h1computesym(N)
    % Comptues the wavelet coefficients of the orthonormal wavelet with N
    % vanishing moments and close to linear phase. This makes the wavelet
    % almost symmetric. These wavelets are called 'symlets'
    %
    % This function relies on matlabs wavelet coefficients. In the next version
    % this will be changed 
    currDWTmode = dwtmode('status', 'nodisp');
    dwtmode('per','nodisp');
    nu = 7;
    n = 2^nu;
    x = zeros([1,n]);
    x(ceil(N/2)) = 1;
    
    S = [2^(nu-1); 2^(nu-1); n]; % compute the S given by wavedec
    wave_name = sprintf('sym%d', N);
    
    y = waverec(x, S, wave_name);
    if (mod(N,2) == 1) % is odd
        g0 = y(1:2*N);
    else % is even 
        g0 = [y(end), y(1:2*N-1)];
    end
    
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); 
    h1=fliplr(g1);
    
    dwtmode(currDWTmode, 'nodisp');

end
