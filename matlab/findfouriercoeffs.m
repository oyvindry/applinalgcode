function [ank,bnk] = findfouriercoeffs(n, k, T)
    ank=0; bnk=0;
    if k > 0
        [ankprev,bnkprev] = findfouriercoeffs(n, k-1, T)
        ank = -k*T*bnkprev/(2*pi*n);
        bnk = -T^k/(pi*n) + k*T*ankprev/(2*pi*n);
    end