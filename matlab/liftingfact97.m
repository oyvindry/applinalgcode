function [lambdas, alpha, beta]=liftingfact97()
    [h0, h1] = h0h1compute97() % Should have 9 and 7 filter coefficients.
    h00 = h0(1:2:9);
    h01 = h0(2:2:9);
    h10 = h1(1:2:7);
    h11 = h1(2:2:7); 
    
    lambdas=zeros(1,4);
    
    lambdas(1) = -h00(1)/h10(1);
    h00(1:5) = h00(1:5) + conv(h10(1:4),[lambdas(1) lambdas(1)]);
    h01(1:4) = h01(1:4) + conv(h11(1:3),[lambdas(1) lambdas(1)]); 
    
    lambdas(2) = -h10(1)/h00(2);
    h10(1:4) = h10(1:4)+conv(h00(2:4),[lambdas(2) lambdas(2)]);
    h11(1:3) = h11(1:3)+conv(h01(2:3),[lambdas(2) lambdas(2)]);
    
    lambdas(3) = -h00(2)/h10(2);
    h00(2:4) = h00(2:4)+conv(h10(2:3),[lambdas(3) lambdas(3)]);
    h01(2:3) = h01(2:3)+conv(h11(2:2),[lambdas(3) lambdas(3)]); 
    
    lambdas(4) = -h10(2)/h00(3);
    h10(2:3) = h10(2:3)+conv(h00(3:3),[lambdas(4) lambdas(4)]);
    
    alpha = h00(3);
    beta  = h11(2);
end

function [h0,h1]=h0h1compute97()
  N=4;
  vals=computeQN(N);
   
  rts=roots(vals)';
  rts1=rts(find(abs(imag(rts))>0.001)); % imaginary roots
  rts2=rts(find(abs(imag(rts))<0.001)); % real roots
  
  h0=1;
  for rt=rts1
    h0=conv(h0,[-rt 1]);
  end
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  
  g0=1;
  for rt=rts2
    g0=conv(g0,[-rt 1]);
  end
  for k=1:(N/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  
  h0=real(h0);
  g0=real(g0);	
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end
