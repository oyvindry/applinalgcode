x=dwt_impl((1:63)', 'spline4.4',  1, 'bd', 'bd_pre', 0, 0, 'resolution'); x'
x=idwt_impl(x, 'spline4.4',  1, 'bd', 'bd_pre', 0, 0, 'resolution');
x'

x = double(imread('images/lena.png'));
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


x=dwt_impl((1:64)', 'db2',  1, 'bd', 'bd_pre', 0, 0, 'resolution');
x'
x=idwt_impl(x, 'db2', 1,  'bd', 'bd_pre', 0, 0, 'resolution');
x'


x=dwt_impl(eye(10), 'db2', 1, 'bd', 'none', 0, 0, 'time');
x
    






def filter_impl(t, x, bd_mode):
    szx = shape(x)
    N = szx[0]
    n = int(prod(szx[1:]))
    y =  reshape(x, (N, n))
    tlen = len(t); N0 = int((tlen - 1)/2)
    w = 0
    
    if bd_mode.lower() == 'symm':
        w = concatenate([ y[N0:0:(-1)], y, y[(N-2):(N - N0 - 2):(-1)] ])
    elif bd_mode.lower() == 'per':
        w = concatenate([ y[(N - N0):], y, y[:N0]])
    elif bd_mode.lower() == 'none' or bd_mode.lower() == 'bd':
        w = concatenate( (zeros(N0, n), y, zeros(N0, n)) )
    for k in range(n):
        z = convolve(t, w[:, k])
        y[:, k] = z[(2*N0):(len(z)-2*N0)]
    x[:] = reshape(y, szx)