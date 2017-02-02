function X=forw_comp_rev_DFT2(f, invf, threshold)
    X = CreateExcerpt();
    M = size(X, 1);
    N = size(X, 2);
    
    X = tensor_impl(X, f, f);
    tot = prod(size(X));

    thresholdmatr = (abs(X) >= threshold);
    zeroedout = tot - sum(sum(sum(thresholdmatr)));
    X = X.*thresholdmatr;
    X = tensor_impl(X, invf, invf);
    X = abs(X);
    X = mapto01(X);
    X = X*255;
    disp sprintf('%f percent of samples zeroed out', 100*zeroedout/tot);
    