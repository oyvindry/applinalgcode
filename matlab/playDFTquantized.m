function playDFTquantized(n)
    [x, fs] = audioread('sounds/castanets.wav');
    N = 32;
    sz = size(x, 1);
    left = sz;
    numzeroed = 0;
    while left >= N
        next = sz + 1 - left;
        y = fft(x(next:(next + N - 1), :));
        y = y/2^n;
        y = round(y);
        y = y*2^n;
        x(next:(next + N - 1), :) = ifft(y);
        left = left - N;
    end
    x = real(x)/max(max(abs(x)));
    playerobj = audioplayer(x, fs);
    playblocking(playerobj);