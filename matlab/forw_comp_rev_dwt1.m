function [x, fs]=forw_comp_rev_dwt1(m, wave_name, lr)
    lowres = 1;
    if nargin >= 3
        lowres = lr;
    end
    [x, fs] = audioread('sounds/castanets.wav');
    N = 2^17;
    x = x(1:N,1);
    x = dwt_impl(x, wave_name, m);
    if lowres==1
        x((N/2^m+1):end) = 0;
    else
        x(1:(N/2^m)) = 0;
    end
    x = idwt_impl(x, wave_name, m);
    x = x/max(max(abs(x)));
end