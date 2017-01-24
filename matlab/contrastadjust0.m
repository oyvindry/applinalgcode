function Z=contrastadjust0(X,n)
  Z = X/255; % Maps the pixel values to [0,1]
  Z = atan(n*(Z-1/2))/(2*atan(n/2)) + 1/2;
  Z = Z*255; % Maps the values back to [0,255]