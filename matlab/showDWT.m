function showDWT(m, wave_name, invf, lowres)
    img = double(imread('images/lena.png', 'png'));
    img = DWT2Impl(img, m, wave_name);
    if lowres
        M =size(img, 1); N = size(img, 2);
        tokeep=img(1:(M/(2^m)), 1:(N/(2^m)), :);
        img=zeros(size(img));
        img(1:(M/(2^m)), 1:(N/(2^m)), :)=tokeep;
    else
        sz = size(img);
        sz(1) = sz(1)/2^m; sz(2) = sz(2)/2^m;
        img(1:sz(1), 1:sz(2), :) = 0;
    end
    img = IDWT2Impl(img, m, invf);
    imshow(uint8(255*mapto01(img)));
  
  
  
  
