function img=forw_comp_rev_dwt2(img, m, wave_name, lr)
    lowres = 1;
    if nargin >= 4
        lowres = lr;
    end
    M = size(img, 1);
    N = size(img, 2);
    img = dwt_impl(img, wave_name, m);
    if lowres==1
        tokeep = img(1:(M/(2^m)), 1:(N/(2^m)), :);
        img=zeros(size(img));
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = tokeep;
    else
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = 0;
    end
    img = idwt_impl(img, wave_name, m);
    img = map_to_01(img);
    img = img*255;
end