function playDFTthreshold(threshold)
    [x, fs] = audioread('sounds/castanets.wav');
    N = 32;
    sz = size(x, 1);
    left = sz;
    numzeroed = 0;
    while left >= N
        next = sz + 1 - left;
        y = fft(x(next:(next + N - 1), :)); 
        numzeroed = numzeroed + sum(sum((abs(y)<threshold)));
        y = (abs(y)>=threshold).*y;
        x(next:(next + N - 1), :) = ifft(y);
        left = left - N;
    end
    100*numzeroed/prod(size(x))
    x = real(x)/max(max(abs(x)));
    playerobj = audioplayer(x, fs);
    playblocking(playerobj);