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

