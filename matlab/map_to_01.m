function Z=map_to_01(X)
  minval = min(min(min(X)));
  maxval = max(max(max(X)));
  Z = (X - minval)/(maxval-minval);
end