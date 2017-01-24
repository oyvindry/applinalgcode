function y = FFTImpl(x, FFTKernel, forward)
    fwd = 1;
    if nargin >= 3
        fwd = forward;
    end
    x = bitreverse(x);
    [N, n] = size(x);
    y = zeros(N, n);
    for s2 = 1:size(x, 2)
        y(:, s2) = FFTKernel(x(:,s2), fwd);
    end
    if ~fwd
        y = y/N;
    end