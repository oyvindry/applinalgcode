x = double(imread('images/lena.png'));
x=dwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre'); 
x=idwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre');
imshow(uint8(x))



    