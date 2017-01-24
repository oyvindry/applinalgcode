function x=bitreverse(x)
    % bitreverse(x) returns x in bit-reversed order
    N = size(x, 1);
    j=0;
    for i = 0:2:(N/2 - 2)
        if (j > i)         
            temp = x(j+1, :); x(j+1, :) = x(i+1, :); x(i+1, :) = temp;
            temp = x(N-i, :); x(N-i, :) = x(N-j, :); x(N-j, :) = temp;
        end
        temp = x(i+2, :); x(i+2, :) = x(j + N/2 + 1, :); x(j + N/2 + 1, :) = temp;
        m = N/4;
        while (m >= 1 && j >= m) 
            j = j - m;
            m = m/2;
        end
        j = j + m;
    end