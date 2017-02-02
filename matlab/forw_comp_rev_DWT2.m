def img=forw_comp_rev_DWT2(m, wave_name, lr)
    lowres = 1;
    if nargin >= 3
        lowres = lr;
    end
    img = CreateExcerpt();
    M = size(img, 1);
    N = size(img, 2);
    img = DWT2Impl(img, m, wave_name);
    if lowres==1
        tokeep = img(1:(M/(2^m)), 1:(N/(2^m)), :);
        img=zeros(size(img));
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = tokeep;
    else
        img(1:(M/(2^m)), 1:(N/(2^m)), :) = 0;
    end
    img = IDWT2Impl(img, m, wave_name);
    img = mapto01(img);
    img = img*255;