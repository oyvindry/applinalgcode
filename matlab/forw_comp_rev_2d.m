function X=forw_comp_rev_2d(X, f, invf, threshold)
    X = tensor2_impl(X, f, f, 'symm');
    tot = prod(size(X));
  
    thresholdmatr = (abs(X) >= threshold);
    zeroedout = tot - sum(sum(sum(thresholdmatr)));
    X = X.*thresholdmatr;
    X = tensor2_impl(X, invf, invf, 'symm');
    X = real(X);
    X = map_to_01(X);
    X = X*255;
    disp(sprintf('%f percent of samples zeroed out', 100*zeroedout/tot));
end