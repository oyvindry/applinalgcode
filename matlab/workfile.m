L = 4;
N = 2^L;
nres = 2;
vm = 2;
mode = 'bd';
dual = 0;


wave_name = sprintf('db%d', vm);
wave_name_x = sprintf('db%dx', vm);

filter = liftingfactortho(2, 0); % 0 is 'db', 1 is 'sym'
filter.A_L_pre;



%x = rand([N,1]);
x = ones([N,1]);

%y = DWTImpl(x, 1, wave_name, mode, dual);
z = FWT_CDJV(x, 2, 2, 'myMode');



%y = IDWTImpl(DWTImpl(x,nres, wave_name,   mode, dual), nres, wave_name,   mode, dual);
%z = IDWTImpl(DWTImpl(x,nres, wave_name_x, mode, dual), nres, wave_name_x, mode, dual);
%
%norm(x-y)
%norm(x-z)




