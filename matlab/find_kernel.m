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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            Find Kernels                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           Kernel functions                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
function x=dwt_kernel_pwl0(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
end

function x=dwt_kernel_pwl0_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
end    
    
function x=idwt_kernel_pwl0(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_odd_symm(0.5, x, bd_mode);
end

function x=idwt_kernel_pwl0_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_even_symm(-0.5, x, bd_mode);
end
    
function x=dwt_kernel_pwl2(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
    x = lifting_even_symm(0.25, x, bd_mode);
end

function x=dwt_kernel_pwl2_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
    x = lifting_odd_symm(-0.25, x, bd_mode);
end
    
function x=idwt_kernel_pwl2(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(-0.25, x, bd_mode);
    x = lifting_odd_symm(0.5, x, bd_mode);
end

function x=idwt_kernel_pwl2_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(0.25, x, bd_mode);
    x = lifting_even_symm(-0.5, x, bd_mode);
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
    
function x=dwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
    if strcmpi(bd_mode, 'bd_pre')
        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
    end
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        y1 = filters.AL'*x(1:size(filters.AL,1));
        y2 = filters.AR'*x((N-size(filters.AR,1)+1):N);
    
        x(1:2:N, :) = x(1:2:N, :)/filters.alpha;
        x(2:2:N, :)=x(2:2:N, :)/filters.beta;
        for stepnr = size(filters.lambdas,1):(-2):2
            x = lifting_odd(filters.lambdas(stepnr,2), filters.lambdas(stepnr,1), x, bd_mode);
            x = lifting_even(filters.lambdas(stepnr-1,2), filters.lambdas(stepnr-1,1), x, bd_mode);
        end
  
        if stepnr == 3
            x = lifting_odd(filters.lambdas(1,2), filters.lambdas(1,1), x, bd_mode);
        end
    else
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
        x(1:size(filters.AL,2)) = x(1:size(filters.AL,2)) + y1;
        x((N-size(filters.AR,2)+1):N) = x((N-size(filters.AR,2)+1):N) + y2;
    end
end

function x=dwt_kernel_ortho_dual(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
    if strcmpi(bd_mode, 'bd_pre')
        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
    end
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

function x=idwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if ( strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre') )
       y1 = filters.AL*x(1:size(filters.AL,2));
       y2 = filters.AR*x((N-size(filters.AR,2)+1):N);
    end

    stepnr = 1;
    if mod(size(filters.lambdas, 1), 2) == 1
        x = lifting_even(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, bd_mode);
        stepnr = stepnr + 1;
    end
  
    while stepnr < size(filters.lambdas, 1)
        x = lifting_odd(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, bd_mode);
        stepnr = stepnr + 1;
        x = lifting_even(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, bd_mode);
        stepnr = stepnr + 1;
    end
 
    x(1:2:N, :)=x(1:2:N, :)/filters.alpha;
    x(2:2:N, :)=x(2:2:N, :)/filters.beta;

    if ( strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre') )
        x(1:size(filters.AL,1)) = x(1:size(filters.AL,1)) + y1;
        x((N-size(filters.AR,1)+1):N) = x((N-size(filters.AR,1)+1):N) + y2;
        if strcmpi(bd_mode, 'bd_pre')
            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
        end
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
        if strcmpi(bd_mode, 'bd_pre')
            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           Lifting Steps                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
function x=lifting_even(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    if strcmpi(bd_mode, 'per')
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
end
  
function x=lifting_odd(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    x(2:2:(N-1), :) = lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if strcmpi(bd_mode, 'per')
        x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre')
        x(N, :) = x(N, :) + lambda2*x(N-1, :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Support functions non-orthogonal                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Support functions orthogonal                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filter=getDBfilter(vm, type)
    %cur_dir =  fileparts(mfilename('fullpath'));
    %
    %if (type == 0)
    %    filename = sprintf('%s/var/DB%d.mat', cur_dir, vm);
    %else
    %    filename = sprintf('%s/var/sym%d.mat', cur_dir, vm);
    %end
    %
    %if (exist(filename) == 2);
    %    load(filename);
    %else
        filter = liftingfactortho(vm, type);
        %global AL AR;
        %AL = filter.AL;
        %AR = filter.AR
    %    
    %    if (exist(strcat(cur_dir, '/var')) ~= 7)
    %        mkdir(strcat(cur_dir,'/var'));
    %    end
    %    
    %    save(filename, 'filter');
    %end
end

function vals=computeQN(N)
    % Compute the coefficients in Q^(N)((1-cos(w))/2)
    k=0:(N-1);
    QN = 2*factorial(N+k-1)./(factorial(k).*factorial(N-1));
    vals=zeros(1,2*N-1);
    vals=QN(1);
    start=1;
    for k=2:N
        start=conv(start,[-1/4 1/2 -1/4]);
        vals=[0 vals 0]+QN(k)*start;
    end
end

function [h0, h1, g0, g1]=h0h1computeortho(N)
    % Comptues the wavelet coefficients of the orthonormal Daubechies wavelet
    % N vanishing moments and with minimum phase   
    vals=computeQN(N);
    rts=roots(vals)';
    rts1=rts(find(abs(rts)>1));

    g0=1;
    for rt=rts1
        g0=conv(g0,[-rt 1]);
    end
    g0 = real(g0);
    K=sqrt(vals(1)*(-1)^(length(rts1))/abs(prod(rts1)));
    g0=K*g0;
    for k=1:N
        g0=conv(g0,[1/2 1/2]);
    end
    h0=fliplr(g0);
    
    % Ensure that sum(h0) == sqrt(2).
    % This is necessary in order for the wavelet to have integral 1 and not -1. 
    if(sum(h0) < -1) 
        h0 = -h0;
        g0 = -g0;
    end

    g1=h0.*(-1).^(0:(length(g0)-1)); 
    h1=fliplr(g1);

end

function [h0, h1, g0, g1]=h0h1computesym(N)
    % Comptues the wavelet coefficients of the orthonormal wavelet with N
    % vanishing moments and close to linear phase. This makes the wavelet
    % almost symmetric. These wavelets are called 'symlets'
    %
    % This function relies on matlabs wavelet coefficients. In the next version
    % this will be changed 
    currDWTmode = dwtmode('status', 'nodisp');
    dwtmode('per','nodisp');
    nu = 7;
    n = 2^nu;
    x = zeros([1,n]);
    x(ceil(N/2)) = 1;
    
    S = [2^(nu-1); 2^(nu-1); n]; % compute the S given by wavedec
    wave_name = sprintf('sym%d', N);
    
    y = waverec(x, S, wave_name);
    if (mod(N,2) == 1) % is odd
        g0 = y(1:2*N);
    else % is even 
        g0 = [y(end), y(1:2*N-1)];
    end
    
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); 
    h1=fliplr(g1);
    
    dwtmode(currDWTmode, 'nodisp');
end

function filters=liftingfactortho(N, type, debug_mode)
    % Computes the filter coefficients of orthonormal wavelets with N vanishing
    % moments.
    %
    % N:    Number of vanishing moments
    % type: The type of orthonormal wavelet.
    %       0: Daubechies wavelets with minimum phase (default)  
    %       1: Symlets - wavelets with close to linear phase (almost symetric) 
    % debug_mode: Wheter or not this function should be ran in debug mode
    % 
    if (nargin <  3)
        debug_mode = 0;
    end
    if (nargin == 1)
        type = 0;
    end
        
    % First the right edge
    if (type == 0)
        [h0, h1, g0, g1] = h0h1computeortho(N);
    elseif (type == 1)
        [h0, h1, g0, g1] = h0h1computesym(N);
    end
    h0 = flip(h0);
    h1 = flip(h1);
    g0 = flip(g0);
    g1 = flip(g1);
    filters = liftingstepscomputeortho(h0, h1);
     
    [W, A_pre, A_pre_inv] = bw_compute_left(h0, g0, debug_mode); % Lower right (3N-1)x(2N) matrix
    
    WR = zeros(size(W));
    for k=1:N
        WR(:,[2*k-1 2*k]) = W(size(W,1):(-1):1,2*N+1-[2*k 2*k-1]); 
    end
     
    % Then the left edge
    h0 = flip(h0);
    h1 = flip(h1);
    g0 = flip(g0);
    g1 = flip(g1);
     
    filters = liftingstepscomputeortho(h0, h1);
    filters.A_R_pre = fliplr(flipud(A_pre));
    filters.A_R_pre_inv = fliplr(flipud(A_pre_inv));
    [WL,A_pre, A_pre_inv] = bw_compute_left(h0, g0, debug_mode); % Upper left (3N-1)x(2N) matrix
    filters.A_L_pre = A_pre;
    filters.A_L_pre_inv = A_pre_inv;
     
    % Compute the left and right parts of the IDWT for boundary handling
    M = 6*N;
    seg1 = zeros(M); % One bigger than is actually needed
     
    filters.AL = zeros(size(WL));
    filters.AR = zeros(size(WR));
    
    for k=0:(M-1)
        x = zeros(M,1);
        x(k+1) = 1;
        seg1(:,k+1) = idwt_kernel_ortho(x, filters, 'bd');
    end
    
    [w1, w2] = size(WL);
    filters.AL=WL-seg1(1:w1,1:w2);
    filters.AR=WR-seg1((M-w1+1):M,(M-w2+1):M);
end

function [W,A_pre,A_pre_inv]=bw_compute_left(h0, g0, debugMode)
    if (nargin < 3) 
        debugMode = 0;
    end
    N = length(h0)/2;
    % We use the symbolic toolbox for computations with rational numbers 
    C_c = sym(zeros(N));
    C_c(1,1) = sym(1);

    L = 2/sym((2*N-1));
    C_c(2,2) = L;

    for n = 1:(N-2) % fill out column n+2. n is as in paper
        betaval = (sym(n^2)/sym(4*n^2-1))*(1-sym(n^2)/sym((2*N-1)^2)) ;
        C_c(2:N,n+2) = L*C_c(1:(N-1),n+1);
        C_c(:,n+2) = C_c(:,n+2) - betaval*C_c(:,n);
    end

    C_0 = sym(zeros(2*N-1,N));
    for k=0:(N-1)
        C_0(:,k+1) = sym((-N+1):(N-1)).^k;
    end

    C = C_0*C_c;
    m = C'*C;
    for k=1:N
        C(:,k) = C(:,k)/sym(sqrt(m(k,k)));
    end

    % dimension of C is (2*N-1)xN
    C = double(C); % Use numerically from now on
    %C'*C % Check for orthogonality

    % For the left edge we need the ((-N+1):(3N-2))x((-2N+2):(2N-2)) segment of
    % G. The first 2N-1 rows are for X, the remaining 2N-1 rows are for Z.
    % This equals the ((N-1):(5N-4))x(0:(4N-3)) segment of G, due to the
    % Toeplitz structure.
    
    seg1 = zeros(2*N-1);
    seg2 = zeros(2*N-1);
    for i =1:N-1
        seg1(1:2*i, i) = g0((end-2*i+1):end); 
    end
    for i = 0:N-1
        seg1((1+2*i):end, N+i) = g0(1:(end - 2*i - 1));
        seg2(1:(2*i+1), N+i) = g0((end - 2*i):end);
    end
    
    % Theorem 3.3
    X_e = C'*seg1*C; % dimensions NxN
    Z_e =    seg2*C; % dimensions (2N-1)xN
    
    % Step 1: Make phi-supports staggered.
    [Q,R] = qr((flipud(C(N:end,:)))');
    P = fliplr(Q);
    invP = flipud(Q');

    % Theorem 3.4
    X_e = invP*X_e*P;
    Z_e = Z_e*P;
    A_pre_inv = C(N:end,:)*P;

    % Step 2: Orthogonalize phi-functions
    ls = eye(N^2)-kron(X_e',X_e');
    rs = reshape(Z_e'*Z_e, [N^2, 1]); % From matrix to vector
    Y  = reshape(ls\rs, [N, N]);      % From vector to matrix
    P = ortho_from_gramm(Y, N);
    
    % Theorem 3.4
    X_e = inv(P)*X_e*P;
    Z_e = Z_e*P;
    A_pre_inv = A_pre_inv*P;
    A_pre = inv(A_pre_inv);
    
    if (debugMode == 1)
        % Step 2: Orthogonalize phi-functions
        ls = eye(N^2)-kron(X_e',X_e');
        rs = reshape(Z_e'*Z_e, [N^2, 1]); % From matrix to vector
        Y  = reshape(ls\rs, [N, N]);      % From vector to matrix
        P = ortho_from_gramm(Y, N);
        
        % Theorem 3.4
        X_e = inv(P)*X_e*P;
        Z_e = Z_e*P;
        A_pre_inv = A_pre_inv*P;
        A_pre = inv(A_pre_inv);
    end
    
    % Step 3: Project the psi-functions onto the orthogonal complement of the,
    % phi-functions. 
    X_o = eye(N) - X_e*X_e';
    Z_o = -Z_e*X_e';

    % Step 4: Make psi-supports staggered
    [Q,R] = qr((flipud(Z_o(1:2:(2*N-1),:)))');
    P = fliplr(Q);
    
    % Theorem 4.1
    X_o = X_o*P;
    Z_o = Z_o*P;

    % Step 5: Orthogonalize psi-functions
    Y = X_o'*X_o + Z_o'*Z_o; % Gram matrix
    P=ortho_from_gramm(Y, N);
    
    % Theorem 4.1
    X_o = X_o*P;
    Z_o = Z_o*P;
    
    W = zeros(3*N-1,2*N);
    W(:, 1:2:(2*N-1)) = [X_e; Z_e];
    W(:, 2:2:(2*N)) = [X_o; Z_o];
end    

function P=ortho_from_gramm(Y, N)
    xmatr=zeros(N);
    for k=2:N
        xmatr(1:(k-1),k) = Y(1:(k-1),1:(k-1))\Y(1:(k-1),k);
    end
    P = eye(N) - xmatr;             % The g_n are now orthogonal
    P = P*diag(sqrt(diag(P'*Y*P)).^(-1)); % The g_n are now orthonormal
end

function filters = liftingstepscomputeortho(h0, h1)
    %global wavlib_lambdas wavlib_alpha wavlib_beta h0 h1
    
    stepnr=1;
    len1=length(h0)/2; len2=len1;
    filters.lambdas=zeros(len1+1,2);
    if mod(len1,2)==0
        h00=h0(1:2:length(h0));
        h01=h0(2:2:length(h0));
        h10=h1(1:2:length(h1));
        h11=h1(2:2:length(h1));
  
        lambda1=-h00(1)/h10(1);
        h00=h00+lambda1*h10; 
        h01=h01+lambda1*h11;
        start1=2; end1=len1; len1=len1-1; start2=1; end2=len2;
        filters.lambdas(stepnr,:)=[lambda1 0];
    else
        h00=h0(2:2:length(h0));
        h01=h0(1:2:length(h0));
        h10=h1(2:2:length(h1));
        h11=h1(1:2:length(h1));
    
        lambda1=-h10(len1)/h00(len1); 
        h10=h10+lambda1*h00; 
        h11=h11+lambda1*h01;
        start2=1; end2=len2-1; len2=len2-1; start1=1; end1=len1;
        filters.lambdas(stepnr,:)=[0 lambda1];
    end
  
    %[h00 h01; h10 h11]
    %conv(h00,h11)-conv(h10,h01)
    stepnr=stepnr+1;

    %[h00 h01; h10 h11]
    %conv(h00,h11)-conv(h10,h01)
    while len2>0 % Stop when the second element in the first column is zero
        if len1>len2 % Reduce the degree in the first row. 
            lambda1=-h00(start1)/h10(start2);
            lambda2=-h00(end1)/h10(end2);
            h00(start1:end1)=h00(start1:end1)+conv(h10(start2:end2),[lambda1 lambda2]);
            h01(start1:end1)=h01(start1:end1)+conv(h11(start2:end2),[lambda1 lambda2]);
            start1=start1+1; end1=end1-1; len1=len1-2;
        else % reduce the degree in the second row. 
            lambda1=-h10(start2)/h00(start1);
            lambda2=-h10(end2)/h00(end1);
            h10(start2:end2)=h10(start2:end2)+conv(h00(start1:end1),[lambda1 lambda2]);
            h11(start2:end2)=h11(start2:end2)+conv(h01(start1:end1),[lambda1 lambda2]);
            start2=start2+1; end2=end2-1; len2=len2-2;
        end
        filters.lambdas(stepnr,:)=[lambda1 lambda2];
        stepnr=stepnr+1;
    
        %[h00 h01; h10 h11]
        %conv(h00,h11)-conv(h10,h01)
    end
  
    % Add the final lifting, and alpha,beta
    filters.alpha=sum(h00);
    filters.beta=sum(h11);
    lastlift=-sum(h01)/filters.beta;
    if mod(length(h0)/2,2)==0
        filters.lambdas(stepnr,:)=[0 lastlift];
    else
        filters.lambdas(stepnr,:)=[lastlift 0];
    end
    %[h00 h01; h10 h11]
end

