function playDFT(L, lower)
    [x, fs] = audioread('sounds/castanets.wav');
    N=32;
    sz = size(x, 1);
    left = sz;
    while left >= N
        next = sz + 1 - left;
        y = fft(x(next:(next + N - 1), :)); 
        if lower
            y((L+2):(N-L), :) = 0;
        else
            y(1:(N/2 - L), :) = 0;
            y((N/2 + L + 1):N, :) = 0;
        end
        x(next:(next + N - 1), :) = ifft(y);
        left = left - N;
    end
    x = real(x)/max(max(abs(x)));
    playerobj = audioplayer(x, fs);
    playblocking(playerobj);