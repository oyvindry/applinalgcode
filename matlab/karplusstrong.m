function karplusstrong(x_init,f_s)
    p=size(x_init,1)-1;
    antsec = 10;
    numsamples = f_s*antsec;
    z=zeros(numsamples,1);
    z(1:(p+1)) = x_init;
    for k=(p+2):numsamples
        z(k) = 0.5*(z(k-p)+z(k-p-1));
    end
    playerobj = audioplayer(z,f_s);
    playblocking(playerobj)
end