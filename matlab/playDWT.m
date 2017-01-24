function playDWT(m, f, invf, lowres)
    [x fs] = audioread('sounds/castanets.wav'); 
    N=2^17;
    x = DWTImpl(x(1:N,:), m, f);
    if lowres
        x((N/2^m+1):N, :) = 0;
    else
        x(1:(N/2^m), :) = 0;
    end
    x = IDWTImpl(x, m, invf);
    playerobj = audioplayer(x, fs);
    playblocking(playerobj);