function x=recover_x(y,N)
    s= size(y,1)/2;
    A=zeros(s);
    for k=1:s
        A(:,s-k+1) = y(k:(k+s-1));
    end
    b=-y((s+1):(2*s));
    w=A\b;
    v = [1; w; zeros(N-(s+1),1)];
    inds = find(abs(ifft(v)) <= 0.00001);
    x_S = exp(-2*pi*i*((0:(s-1))')*(inds-1)'/N)\y(1:s);
    x = zeros(N,1);
    x(inds)=x_S;