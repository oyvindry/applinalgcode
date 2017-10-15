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
    
    % Ensuring integral is positive - This part of the code requiere some more
    % testing
    if (sum(g0) < 0)
        g0 = -g0;
    end
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); % It seems to me that this should be 
                                     % multiplied by -1 
    h1=fliplr(g1);
end

