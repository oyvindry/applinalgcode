function playDWT(m, wave_name, lowres)
    [x fs] = audioread('sounds/castanets.wav'); 
    N=2^17;
    x = dwt_impl(x(1:N,:), m, wave_name);
    if lowres
        x((N/2^m+1):N, :) = 0;
    else
        x(1:(N/2^m), :) = 0;
    end
    x = dwt_impl(x, m, wave_name, 0);
    playerobj = audioplayer(x, fs);
    playblocking(playerobj);