function x = mp3_reverse_fbt(z)
    Ns = length(z)/32;
    x = zeros(32*Ns, 1);
    D = mp3_dtable();          % The synthesis window.
    V = zeros(1024,1);
    % The 64x32 matrix N
    N = cos((16+((0:63)'))*(2*(0:31)+1)*pi/64);
  
    U = zeros(512,1);
    for n = 1:Ns;
        V(65:1024) = V(1:(1024-64));
        V(1:64) = N*z((1+(n-1)*32):(n*32));
        for i = 0:7
            U((i*64 + 1):(i*64 + 32))   = V((i*128 + 1):(i*128 + 32));
            U((i*64 + 33):((i + 1)*64)) = V((i*128 + 97):((i+1)*128));
        end
        W = U.*D;
        for i = 0:15
            x(((n-1)*32 + 1):(n*32)) = x(((n-1)*32 + 1):(n*32)) ...
                                       + W((32*i + 1):(32*(i + 1)));
        end 
    end