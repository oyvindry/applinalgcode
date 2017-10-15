function f=find_kernel(wave_name, forward, dual, transpose)
    if transpose
        forward = ~forward;
        dual = ~dual;
    end
    if forward
        if dual
            f = find_kernel_dwt_dual(wave_name);
        else
            f = find_kernel_dwt(wave_name);
        end
    else
        if dual
            f = find_kernel_idwt_dual(wave_name);
        else
            f = find_kernel_idwt(wave_name);
        end
    end
end

function f= find_kernel_dwt_dual(wave_name)
    if strcmpi(wave_name, 'cdf97')
        f = @dwt_kernel_97_dual;
    elseif strcmpi(wave_name, 'cdf53')
        f = @dwt_kernel_53_dual;
    elseif strcmpi(wave_name, 'pwl0')
        f = @dwt_kernel_pwl0_dual;
    elseif strcmpi(wave_name, 'pwl2')
        f = @dwt_kernel_pwl2_dual;
    elseif strcmpi(wave_name, 'Haar')
        f = @dwt_kernel_haar;
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        if (vm == 1) 
            f = @dwt_kernel_haar;
        else
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        if (vm == 1) 
            f = @dwt_kernel_haar;
        else
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
    elseif strcmpi(wave_name(1:6), 'spline')
        N1 = str2double(wave_name(7));
        N2 = str2double(wave_name(9));
        [h0,h1,g0,g1]=compute_spline_filters(N1, N2);
        f = @(x, bd_mode) dwt_kernel_filters_dual(h0, h1, g0, g1, x, bd_mode);
    end
end
    
function f= find_kernel_dwt(wave_name)    
    if strcmpi(wave_name, 'cdf97')
        % [lambdas, alpha, beta]=liftingfact97()
        f = @dwt_kernel_97;
    elseif strcmpi(wave_name, 'cdf53')
        f = @dwt_kernel_53;
    elseif strcmpi(wave_name, 'pwl0')
        f = @dwt_kernel_pwl0;
    elseif strcmpi(wave_name, 'pwl2')
        f = @dwt_kernel_pwl2;
    elseif strcmpi(wave_name, 'Haar')
        f = @dwt_kernel_haar;
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        if (vm == 1) 
            f = @dwt_kernel_haar;
        else
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        if (vm == 1) 
            f = @dwt_kernel_haar;
        else
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
    elseif strcmpi(wave_name(1:6), 'spline')
        N1 = str2double(wave_name(7));
        N2 = str2double(wave_name(9));
        [h0,h1,g0,g1]=compute_spline_filters(N1, N2);
        f = @(x, bd_mode) dwt_kernel_filters(h0, h1, g0, g1, x, bd_mode);
    end
end

function f = find_kernel_idwt_dual(wave_name)
    if strcmpi(wave_name, 'cdf97')
        f = @idwt_kernel_97_dual;
    elseif strcmpi(wave_name, 'cdf53')
        f = @idwt_kernel_53_dual;
    elseif strcmpi(wave_name, 'pwl0')
        f = @idwt_kernel_pwl0_dual;
    elseif strcmpi(wave_name, 'pwl2')
        f = @idwt_kernel_pwl2_dual;
    elseif strcmpi(wave_name, 'Haar')
        f = @idwt_kernel_haar;
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        if (vm == 1)
            f = @idwt_kernel_haar;
        else
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        if (vm == 1)
            f = @idwt_kernel_haar;
        else
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
    elseif strcmpi(wave_name(1:6), 'spline')
        N1 = str2double(wave_name(7));
        N2 = str2double(wave_name(9));
        [h0,h1,g0,g1]=compute_spline_filters(N1, N2);
        f = @(x, bd_mode) idwt_kernel_filters_dual(h0, h1, g0, g1, x, bd_mode);
    end
end

function f = find_kernel_idwt(wave_name)
    if strcmpi(wave_name, 'cdf97')
        f = @idwt_kernel_97;
    elseif strcmpi(wave_name, 'cdf53')
        f = @idwt_kernel_53;
    elseif strcmpi(wave_name, 'pwl0')
        f = @idwt_kernel_pwl0;
    elseif strcmpi(wave_name, 'pwl2')
        f = @idwt_kernel_pwl2;
    elseif strcmpi(wave_name, 'Haar')
        f = @idwt_kernel_haar;
    elseif (strcmpi(wave_name(1:2), 'db') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end));
        if (vm == 1)
            f = @idwt_kernel_haar;
        else
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(3:end-1));
        filters = liftingfactortho(vm, 0, 1);
        f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
    elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end));
        if (vm == 1)
            f = @idwt_kernel_haar;
        else
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        end
    elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
        vm = str2double(wave_name(4:end-1));
        filters = liftingfactortho(vm, 1, 1);
        f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
    elseif strcmpi(wave_name(1:6), 'spline')
        N1 = str2double(wave_name(7));
        N2 = str2double(wave_name(9));
        [h0,h1,g0,g1]=compute_spline_filters(N1, N2);
        f = @(x, bd_mode) idwt_kernel_filters(h0, h1, g0, g1, x, bd_mode);
    end
end



function [h0,h1,g0,g1]=compute_spline_filters(N1, N2)
  N=(N1+N2)/2;
  vals=computeQN(N);
  
  h0 = 1;
  for k=1:(N1/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  h0 = conv(h0, vals);
  
  g0=1;
  for k=1:(N2/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end

function x=dwt_kernel_filters(H0, H1, G0, G1, x, bd_mode)
    symm = strcmpi(bd_mode, 'symm');
    f0 = H0; f1 = H1;
    x0 = filterS(f0, x, symm);
    x1 = filterS(f1, x, symm);
    x(1:2:end) = x0(1:2:end);
    x(2:2:end) = x1(2:2:end);
end
    
function x=dwt_kernel_filters_dual(H0, H1, G0, G1, x, bd_mode)
    symm = strcmpi(bd_mode, 'symm');
    f0 = G0; f1 = G1;
    x0 = filterS(f0, x, symm);
    x1 = filterS(f1, x, symm);
    x(1:2:end) = x0(1:2:end);
    x(2:2:end) = x1(2:2:end);
end
    
function x=idwt_kernel_filters(H0, H1, G0, G1, x, bd_mode)
    symm = strcmpi(bd_mode, 'symm');
    f0 = G0; f1 = G1;
    x0 = x; x0(2:2:end) = 0;
    x1 = x; x1(1:2:end) = 0;
    x0 = filterS(f0, x0, symm);
    x1 = filterS(f1, x1, symm);
    x = x0 + x1;
end
    
function x=idwt_kernel_filters_dual(H0, H1, G0, G1, x, bd_mode)
    symm = strcmpi(bd_mode, 'symm');
    f0 = H0; f1 = H1;
    x0 = x; x0(2:2:end) = 0;
    x1 = x; x1(1:2:end) = 0;
    x0 = filterS(f0, x0, symm);
    x1 = filterS(f1, x1, symm);
    x = x0 + x1;
end
    
function [lambdas, alpha, beta]=liftingfact97()
    [h0, h1] = h0h1compute97() % Should have 9 and 7 filter coefficients.
    h00 = h0(1:2:9);
    h01 = h0(2:2:9);
    h10 = h1(1:2:7);
    h11 = h1(2:2:7); 
    
    lambdas=zeros(1,4);
    
    lambdas(1) = -h00(1)/h10(1);
    h00(1:5) = h00(1:5) + conv(h10(1:4),[lambdas(1) lambdas(1)]);
    h01(1:4) = h01(1:4) + conv(h11(1:3),[lambdas(1) lambdas(1)]); 
    
    lambdas(2) = -h10(1)/h00(2);
    h10(1:4) = h10(1:4)+conv(h00(2:4),[lambdas(2) lambdas(2)]);
    h11(1:3) = h11(1:3)+conv(h01(2:3),[lambdas(2) lambdas(2)]);
    
    lambdas(3) = -h00(2)/h10(2);
    h00(2:4) = h00(2:4)+conv(h10(2:3),[lambdas(3) lambdas(3)]);
    h01(2:3) = h01(2:3)+conv(h11(2:2),[lambdas(3) lambdas(3)]); 
    
    lambdas(4) = -h10(2)/h00(3);
    h10(2:3) = h10(2:3)+conv(h00(3:3),[lambdas(4) lambdas(4)]);
    
    alpha = h00(3);
    beta  = h11(2);
end

function [h0,h1]=h0h1compute97()
  N=4;
  vals=computeQN(N);
   
  rts=roots(vals)';
  rts1=rts(find(abs(imag(rts))>0.001)); % imaginary roots
  rts2=rts(find(abs(imag(rts))<0.001)); % real roots
  
  h0=1;
  for rt=rts1
    h0=conv(h0,[-rt 1]);
  end
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  
  g0=1;
  for rt=rts2
    g0=conv(g0,[-rt 1]);
  end
  for k=1:(N/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  
  h0=real(h0);
  g0=real(g0);  
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end


function x=dwt_kernel_97_dual(x,bd_mode)
    lambda1 = -0.586134342059950;
    lambda2 = -0.668067171029734;
    lambda3 = 0.070018009414994;
    lambda4 = 1.200171016244178;
    alpha = -1.149604398860250;
    beta = -0.869864451624777;
  
    x(1:2:end, :) =x(1:2:end, :)/alpha;
    x(2:2:end, :)=x(2:2:end, :)/beta;
    x = lifting_even_symm(lambda4, x, bd_mode);
    x = lifting_odd_symm(lambda3, x, bd_mode);
    x = lifting_even_symm(lambda2, x, bd_mode);
    x = lifting_odd_symm(lambda1, x, bd_mode);
end

function x=dwt_kernel_97(x,bd_mode)
    lambda1 = -0.586134342059950;
    lambda2 = -0.668067171029734;
    lambda3 = 0.070018009414994;
    lambda4 = 1.200171016244178;
    alpha = -1.149604398860250;
    beta = -0.869864451624777;
  
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
    x = lifting_odd_symm(-lambda4, x, bd_mode);
    x = lifting_even_symm(-lambda3, x, bd_mode);
    x = lifting_odd_symm(-lambda2, x, bd_mode);
    x = lifting_even_symm(-lambda1, x, bd_mode);
end
    
function x=idwt_kernel_97_dual(x, bd_mode)
    lambda1 = -0.586134342059950;
    lambda2 = -0.668067171029734;
    lambda3 = 0.070018009414994;
    lambda4 = 1.200171016244178;
    alpha = -1.149604398860250;
    beta = -0.869864451624777;
  
    x = lifting_odd_symm(-lambda1, x, bd_mode);
    x = lifting_even_symm(-lambda2, x, bd_mode);
    x = lifting_odd_symm(-lambda3, x, bd_mode);
    x = lifting_even_symm(-lambda4, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
end
    
function x=idwt_kernel_97(x, bd_mode)
    lambda1 = -0.586134342059950;
    lambda2 = -0.668067171029734;
    lambda3 = 0.070018009414994;
    lambda4 = 1.200171016244178;
    alpha = -1.149604398860250;
    beta = -0.869864451624777;
  
    x = lifting_even_symm(lambda1, x, bd_mode);
    x = lifting_odd_symm(lambda2, x, bd_mode);
    x = lifting_even_symm(lambda3, x, bd_mode);
    x = lifting_odd_symm(lambda4, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)/alpha;
    x(2:2:end, :) = x(2:2:end, :)/beta;
end
    
function x=dwt_kernel_53_dual(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x(1:2:end, :) = x(1:2:end, :)/alpha;
    x(2:2:end, :) = x(2:2:end, :)/beta;
    x = lifting_even_symm(lambda2, x, bd_mode);
    x = lifting_odd_symm(lambda1, x, bd_mode);
end
    
function x=dwt_kernel_53(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
    x = lifting_odd_symm(-lambda2, x, bd_mode);
    x = lifting_even_symm(-lambda1, x, bd_mode);
end
    
function x=idwt_kernel_53_dual(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x = lifting_odd_symm(-lambda1, x, bd_mode);
    x = lifting_even_symm(-lambda2, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
end
   
function x=idwt_kernel_53(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x = lifting_even_symm(lambda1, x, bd_mode);
    x = lifting_odd_symm(lambda2, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)/alpha;
    x(2:2:end, :) = x(2:2:end, :)/beta;
end
    
function x=dwt_kernel_pwl0_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
end
    
function x=dwt_kernel_pwl0(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
end
    
function x=idwt_kernel_pwl0_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_even_symm(-0.5, x, bd_mode);
end
    
function x=idwt_kernel_pwl0(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_odd_symm(0.5, x, bd_mode);
end

function x=dwt_kernel_pwl2_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
    x = lifting_odd_symm(-0.25, x, bd_mode);
end
    
function x=dwt_kernel_pwl2(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
    x = lifting_even_symm(0.25, x, bd_mode);
end
    
function x=idwt_kernel_pwl2_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(0.25, x, bd_mode);
    x = lifting_even_symm(-0.5, x, bd_mode);
end

function x=idwt_kernel_pwl2(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(-0.25, x, bd_mode);
    x = lifting_odd_symm(0.5, x, bd_mode);
end
    
function x=dwt_kernel_haar(x, bd_mode)
    x = x/sqrt(2);
    N = size(x, 1);
    if mod(N,2) == 1
        x(1:2, :) = [x(1, :) + x(2, :) - x(N, :); x(1, :) - x(2, :) - x(N, :)];
        x(N, :) = 2*x(N, :);
    else
        x(1:2, :) = [x(1, :) + x(2, :); x(1, :) - x(2, :)];
    end
    for k = 3:2:(N-1)
        x(k:(k+1), :) = [x(k, :) + x(k+1, :); x(k, :) - x(k+1, :)];
    end
end
    
function x=idwt_kernel_haar(x, bd_mode)
    x = x/sqrt(2);
    N = size(x, 1);
    if mod(N,2) == 1
        x(1:2, :) = [x(1, :) + x(2, :) + x(N, :); x(1, :) - x(2, :)];
    else
        x(1:2, :) = [x(1, :) + x(2, :); x(1, :) - x(2, :)];
    end
    for k = 3:2:(N-1)
        x(k:(k+1), :) = [x(k, :) + x(k+1, :); x(k, :) - x(k+1, :)];
    end  
end
    
%% DUAL function
function x=dwt_kernel_ortho_dual(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
%    if strcmpi(bd_mode, 'bd_pre')
%        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
%        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
%    end
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        y1 = filters.AL'*x(1:size(filters.AL,1));
        y2 = filters.AR'*x((N-size(filters.AR,1)+1):N);
        
        x(1:2:N, :) = x(1:2:N, :)*filters.alpha;
        x(2:2:N, :) = x(2:2:N, :)*filters.beta;
        for stepnr = size(filters.lambdas,1):(-2):2
            x = lifting_even(-filters.lambdas(stepnr,1), -filters.lambdas(stepnr,2), x, bd_mode);
            x = lifting_odd(-filters.lambdas(stepnr-1,1), -filters.lambdas(stepnr-1,2), x, bd_mode);
        end
  
        if stepnr == 3
            x = lifting_even(-filters.lambdas(1,1), -filters.lambdas(1,2), x, bd_mode);
        end
    else
        x(1:2:N, :) = x(1:2:N, :)/filters.alpha;
        x(2:2:N, :)=x(2:2:N, :)/filters.beta;
        for stepnr = size(filters.lambdas,1):(-2):2
            x = lifting_odd(filters.lambdas(stepnr,2), filters.lambdas(stepnr,1), x, bd_mode);
            x = lifting_even(filters.lambdas(stepnr-1,2), filters.lambdas(stepnr-1,1), x, bd_mode);
        end
  
        if stepnr == 3
            x = lifting_odd(filters.lambdas(1,2), filters.lambdas(1,1), x, bd_mode);
        end
    end
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(1:size(filters.AL,2)) = x(1:size(filters.AL,2)) + y1;
        x((N-size(filters.AR,2)+1):N) = x((N-size(filters.AR,2)+1):N) + y2;
    end
end





function x=dwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
%    if strcmpi(bd_mode, 'bd_pre')
%        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
%        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
%        %x(1:size(filters.A_L_pre,1))
%    end
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        % Store the boundary coefficients
        y1 = filters.AL'*x(1:size(filters.AL,1));
        y2 = filters.AR'*x((N-size(filters.AR,1)+1):N);
        
        % Why is we dividing on filters.alpha in this step, while multipling on
        % if we do not use boundary wavelets. 
        x(1:2:N, :) = x(1:2:N, :)/filters.alpha;
        x(2:2:N, :) = x(2:2:N, :)/filters.beta;
        for stepnr = size(filters.lambdas,1):(-2):2
            x = lifting_odd(filters.lambdas(stepnr,2), filters.lambdas(stepnr,1), x, bd_mode);
            x = lifting_even(filters.lambdas(stepnr-1,2), filters.lambdas(stepnr-1,1), x, bd_mode);
        end
  
        if stepnr == 3
            x = lifting_odd(filters.lambdas(1,2), filters.lambdas(1,1), x, bd_mode);
        end
    else
        
        % Perform the wavelet transform 
        x(1:2:N, :) = x(1:2:N, :)*filters.alpha;
        x(2:2:N, :) = x(2:2:N, :)*filters.beta;
        for stepnr = size(filters.lambdas,1):(-2):2
            x = lifting_even(-filters.lambdas(stepnr,1), -filters.lambdas(stepnr,2), x, bd_mode);
            x = lifting_odd(-filters.lambdas(stepnr-1,1), -filters.lambdas(stepnr-1,2), x, bd_mode);
        end
  
        if stepnr == 3
            x = lifting_even(-filters.lambdas(1,1), -filters.lambdas(1,2), x, bd_mode);
        end
    end
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        % Add the boundary coefficents 
        x(1:size(filters.AL,2)) = x(1:size(filters.AL,2)) + y1;
        x((N-size(filters.AR,2)+1):N) = x((N-size(filters.AR,2)+1):N) + y2;
    end
end




    
function x=idwt_kernel_ortho_dual(x, filters, bd_mode)
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
       y1 = filters.AL*x(1:size(filters.AL,2));
       y2 = filters.AR*x((N-size(filters.AR,2)+1):N);
    end
    stepnr = 1;
    if mod(size(filters.lambdas, 1), 2) == 1
        x = lifting_odd(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, bd_mode);
        stepnr = stepnr + 1;
    end
  
    while stepnr < size(filters.lambdas, 1)
        x = lifting_even(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, bd_mode);
        stepnr = stepnr + 1;
        x = lifting_odd(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, bd_mode);
        stepnr = stepnr + 1;
    end
 
    x(1:2:N, :) = x(1:2:N, :)*filters.alpha;
    x(2:2:N, :) = x(2:2:N, :)*filters.beta;
    
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(1:size(filters.AL,1)) = x(1:size(filters.AL,1)) + y1;
        x((N-size(filters.AR,1)+1):N) = x((N-size(filters.AR,1)+1):N) + y2;
%        if strcmpi(bd_mode, 'bd_pre')
%            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
%            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
%        end
    end
end


    
function x=lifting_even_symm(lambda, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        assert(mod(N,2) == 0)
    end
    if strcmpi(bd_mode, 'symm')
        x(1, :) = x(1, :) + 2*lambda*x(2, :); % Symmetric extension
    else
        x(1, :) = lambda*(x(2, :) + x(N, :)) + x(1, :);
    end
    x(3:2:(N-1), :) = x(3:2:(N-1), :) + lambda*(x(2:2:(N-2), :) + x(4:2:N, :)); % This saves one multiplication
    if mod(N,2) == 1 % last is odd
        x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
    end
end
    
function x=lifting_odd_symm(lambda, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        assert(mod(N,2) == 0)
    end
    x(2:2:(N-1), :) = x(2:2:(N-1), :) + lambda*(x(1:2:(N-2), :) + x(3:2:N, :)); % This saves one multiplication
    if mod(N,2)==0 % last is even
        if strcmpi(bd_mode, 'symm')
            x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
        else
            x(N, :) = lambda*(x(1, :) + x(N-1, :)) + x(N, :);
        end
    end
end
    







