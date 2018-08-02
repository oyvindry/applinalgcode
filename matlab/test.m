x=DWTImpl((1:63)', 1, 'spline4.4',  'bd', 'bd_pre', 0, 0, 'resolution'); x'
x=IDWTImpl(x, 1, 'spline4.4',  'bd', 'bd_pre', 0, 0, 'resolution');
x'

x = double(imread('../images/lena.png'));
x=dwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre'); 
x=idwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre');
imshow(uint8(x))



        x=dwt_impl((1:64)', 'cdf53',  1, 'bd', 'bd_pre'); 
        x'
        x=idwt_impl(x, 'cdf53', 1, 'bd', 'bd_pre');
        x'
%ok

x=DWTImpl((1:64)', 2, 'pwl2',  'bd', 'bd_pre', 0, 0, 'resolution'); 
x'
x=IDWTImpl(x, 2, 'pwl2',  'bd', 'bd_pre', 0, 0, 'resolution');
x'


x=DWTImpl((1:64)', 1, 'db2',  'bd', 'bd_pre', 0, 0, 'resolution');
x'
x=IDWTImpl(x, 1, 'db2',  'bd', 'bd_pre', 0, 0, 'resolution');
x'


x=DWTImpl(eye(10), 1, 'db2',  'bd', 'none', 0, 0, 'time');
x
    

