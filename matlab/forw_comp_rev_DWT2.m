function img=forw_comp_rev_DWT2(m, wave_name, lr)
    lowres = 1;
    if nargin >= 3
        lowres = lr;
    end
    img = create_excerpt();
    M = size(img, 1);
    N = size(img, 2);
    img = dwt_impl(img, wave_name, n);
    if lowres==1
        tokeep = img(1:(M/(2^m)), 1:(N/(2^m)), :);
        img=zeros(size(img));
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = tokeep;
    else
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = 0;
    end
    img = idwt_impl(img, wave_name, m);
    img = mapto01(img);
    img = img*255;
end