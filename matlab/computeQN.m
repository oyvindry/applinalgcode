function vals=computeQN(N)
  % Compute the coefficients in Q^(N)((1-cos(w))/2)
  k=0:(N-1);
  QN = 2*factorial(N+k-1)./(factorial(k).*factorial(N-1));
  vals=zeros(1,2*N-1);
  vals=QN(1);
  start=1;
  for k=2:N
    start=conv(start,[-1/4 1/2 -1/4]);
    vals=[0 vals 0]+QN(k)*start;
  end