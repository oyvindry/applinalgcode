function X=forw_comp_rev_DCT2(f, invf, threshold)
    X = create_excerpt();
    X = tensor_impl(X, f, f);
    tot = prod(size(X));
  
    thresholdmatr = (abs(X) >= threshold);
    zeroedout = tot - sum(sum(sum(thresholdmatr)));
    X = X.*thresholdmatr;
    X = tensor_impl(X, invf, invf);
    X = mapto01(X);
    X = X*255;
    print sprintf('%f percent of samples zeroed out', 100*zeroedout/tot);
end