L = 7;
N = 2^L;
nres = 2;
vm = 4;
mode = 2;
dual = 0;


wave_name = sprintf('db%d', vm);
wave_name_x = sprintf('db%dx', vm);


x = rand([N,1]);
y = IDWTIm(DWTIm(x,nres, wave_name,   mode, dual), nres, wave_name,   mode, dual);
z = IDWTIm(DWTIm(x,nres, wave_name_x, mode, dual), nres, wave_name_x, mode, dual);

norm(x-y)
norm(x-z)




