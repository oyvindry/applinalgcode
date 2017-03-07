function [x, fs] =forw_comp_rev_DFT(varargin)
    vargs = varargin;
    nargs = length(vargs);
    names = vargs(1:2:nargs);
    values = vargs(2:2:nargs);

    L = 0;
    ind = strmatch('L', names, 'exact');
    if ind>0
        L = values{ind};
    end
    lower = -1;
    ind = strmatch('lower', names, 'exact');
    if ind>0
        lower = values{ind};
    end
    threshold = 0;
    ind = strmatch('threshold', names, 'exact');
    if ind>0
        threshold = values{ind};
    end
    n = 0;
    ind = strmatch('n', names, 'exact');
    if ind>0
        n = values{ind};
    end
    N = 0;
    ind = strmatch('N', names, 'exact');
    if ind>0
        N = values{ind};
    end

    [x, fs] = audioread('sounds/castanets.wav');
    x = x(:,1);
    if N == 0
        N = length(x);
    end
    sz = length(x);
    numzeroed = 0;
    for ind=1:N:sz
        y = fft(x(ind:(ind + N - 1)));
        if lower == 1
            y((L+2):(N-L)) = 0;
        elseif lower == 0
            y(1:(N/2 - L)) = 0;
            y((N/2 + L + 1):N) = 0;
        end
        if threshold ~= 0
            numzeroed = numzeroed + sum((abs(y) < threshold));
            y = y.*(abs(y) >= threshold);
        end
        if n ~= 0
            y = y/2^n;
            y = round(y);
            y = y*2^n;
        end
        x(ind:(ind + N - 1)) = ifft(y);
    end
    x = real(x);
    x = x/max(max(abs(x)));
    if threshold ~= 0
        100*numzeroed/prod(size(x))
    end
end