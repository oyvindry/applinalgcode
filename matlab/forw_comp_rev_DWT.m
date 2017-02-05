function [x,fs]=forw_comp_rev_DWT(m, wave_name, lr)
    lowres = 1;
    if nargin >= 3
        lowres = lr;
    end
    [x, fs] = audioread('sounds/castanets.wav');
    N = 2^17;
    x = x(1:N);
    x = dwt_impl(x, m, wave_name);
    if lowres==1
        x((N/2^m+1):end) = 0;
    else
        x(1:(N/2^m)) = 0;
    end
    x = dwt_impl(x, m, wave_name, 0);
    x = x/max(max(abs(x)));