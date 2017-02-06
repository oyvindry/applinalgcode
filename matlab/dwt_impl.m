function x=dwt_impl(x, nres, wave_name, forward, dim, bd_mode, dual, transpose)
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % nres:      Number of resolutions.
    % wave_name: Name of the wavelet.
    % forward:   Whether to apply the forward or reverse transform. Default: 1
    % dim:       The dimension of the transform (1 for sound, 2 for images). Default: 1
    % bd_mode:   Boundary extension mode. Default: 1
    % dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: 0
    % transpose: Whether the transpose is to be taken. Default: 0
    
    if (~exist('bd_mode','var')) bd_mode = 1; end
    if (~exist('dual','var')) dual  = 0; end
    if (~exist('forward','var')) forward = 1; end
    if (~exist('dim','var')) dim  = 1; end
    if (~exist('transpose','var')) transpose = 0; end
    
    if transpose
        forward = ~forward;
        dual = ~dual;
    end
    if forward
        f = findDWTKernel(wave_name, dual);
        if dim == 2
            x = dwt2_impl(x, nres, f, bd_mode);
        elseif dim == 1
            x = dwt1_impl(x, nres, f, bd_mode);
        end 
    else
        f = findIDWTKernel(wave_name, dual);
        if dim ==2
            x = idwt2_impl(x, nres, f, bd_mode);
        elseif dim == 1
            x = idwt1_impl(x, nres, f, bd_mode);
        end
    end
end    
    
    
    
function x = dwt2_impl(x, nres, f, bd_mode)
    M = size(x, 1); N = size(x, 2); sz = size(x);
    M0 = size(x, 1); N0 = size(x, 2);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = 0:(nres - 1)
        sz2(1) = M; Y2 = zeros(sz2);
        sz1(1) = N; Y1 = zeros(sz1);
        if length(sz1)==1
            Y1=zeros(sz1, 1); Y2=zeros(sz2, 1);
        end
        for n = 1:2^res:N0
            Y2(:, :) = x(1:2^res:M0, n, :);
            x(1:2^res:M0, n, :) = f(Y2, bd_mode);
        end
        for m = 1:2^res:M0
            Y1(:, :) = x(m, 1:2^res:N0, :);
            x(m, 1:2^res:N0, :) = f(Y1, bd_mode);
        end
        M = ceil(M/2); N = ceil(N/2);
    end
    
    x = reorganize_coeffs2(x, nres, 1);   
end    
    
    
function x = dwt1_impl(x, nres, f, bd_mode)
    for res=0:(nres - 1)
        x(1:2^res:end, :) = f(x(1:2^res:end, :), bd_mode);
    end
    x = reorganize_coeffs(x, nres, 1);
end



function x=idwt2_impl(x, nres, f, bd_mode)
    x = reorganize_coeffs2(x, nres, 0);   
    
    M = size(x, 1); N = size(x, 2); sz = size(x);
    sz1 = sz; sz1(1) = [];
    sz2 = sz; sz2(2) = [];
    for res = (nres - 1):(-1):0
        sz1(1) = length(1:2^res:N); sz2(1) = length(1:2^res:M);
        Y1 = zeros(sz1); Y2 = zeros(sz2);
        if length(sz1)==1
            Y1=zeros(sz1, 1); Y2=zeros(sz2, 1);
        end
        for n = 1:2^res:N
            Y2(:, :) = x(1:2^res:M, n, :);
            x(1:2^res:M, n, :) = f(Y2(:, :), bd_mode);
        end
        for m = 1:2^res:M
            Y1(:, :) = x(m, 1:2^res:N, :);
            x(m, 1:2^res:N, :) = f(Y1(:, :), bd_mode);
        end
    end
end    
    
function x = idwt1_impl(x, nres, f, bd_mode)
    x = reorganize_coeffs(x, nres, 0);
    for res = (nres - 1):(-1):0
        x(1:2^res:end, :) = f(x(1:2^res:end, :), bd_mode);
    end
end


function y = reorganize_coeffs(x, nres, forward)
    N = size(x,1);
    y = zeros(size(x));
    inds = 1:2^nres:N;
    lc = length(inds);
    if forward
        y(1:lc, :) = x(inds, :);
    else
        y(inds, :) = x(1:lc, :);
    end
    for res = nres:(-1):1
        inds = (2^(res - 1) + 1):2^res:N;
        lw = length(inds);
        if forward
            y((lc + 1):(lc + lw), :) = x(inds, :);
        else
            y(inds, :) = x((lc + 1):(lc + lw), :);
        end
        lc = lc + lw;
    end
end
    
function Y=reorganize_coeffs2(X, nres, forward)
    M = size(X, 1);
    N = size(X, 2);
    Y = zeros(size(X));
    inds1 = 1:2^nres:M;
    inds2 = 1:2^nres:N;
    lc1 = length(inds1);
    lc2 = length(inds2);
    if forward
        Y(1:lc1, 1:lc2, :) = X(inds1, inds2, :);
    else
        Y(inds1, inds2, :) = X(1:lc1, 1:lc2, :);
    end
    for res = nres:(-1):1
        inds1 = (2^(res - 1) + 1):2^res:M;
        inds2 = (2^(res - 1) + 1):2^res:N;
        lw1 = length(inds1);
        lw2 = length(inds2);
        if forward
            Y((lc1 + 1):(lc1 + lw1), 1:lc2, :) = X(inds1, 1:2^res:M, :);
            Y((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :) = X(inds1, (1+2^(res-1)):2^res:M, :);
            Y(1:lc1, (lc2 + 1):(lc2 + lw2), :) = X(1:2^res:M, inds2, :);
        else
            Y(inds1, 1:2^res:M, :) = X((lc1 + 1):(lc1 + lw1), 1:lc2, :);
            Y(inds1, (1+2^(res-1)):2^res:M, :) = X((lc1 + 1):(lc1 + lw1), (lc2+1):(lc2+lw2), :);
            Y(1:2^res:M, inds2, :) = X(1:lc1, (lc2 + 1):(lc2 + lw2), :);
        end
        lc1 = lc1 + lw1;
        lc2 = lc2 + lw2;
    end    
end

function f= findDWTKernel(wave_name, dual)
    % Find the DWTKernel corresponding to the given wavelet name 
    if dual
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
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(3:end-1));
            filters = liftingfactortho(vm, 0, 1);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end));
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end-1));
            filters = liftingfactortho(vm, 1, 1);
            f = @(x, bd_mode) dwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    else
        if strcmpi(wave_name, 'cdf97')
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
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(3:end-1));
            filters = liftingfactortho(vm, 0, 1);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end));
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end-1));
            filters = liftingfactortho(vm, 1, 1);
            f = @(x, bd_mode) dwt_kernel_ortho(x, filters, bd_mode);
        end
    end
end

function f = findIDWTKernel(wave_name, dual)
    % Find the IDWTKernel corresponding to the given wavelet name
    if dual
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
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(3:end-1));
            filters = liftingfactortho(vm, 0, 1);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end));
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end-1));
            filters = liftingfactortho(vm, 1, 1);
            f = @(x, bd_mode) idwt_kernel_ortho_dual(x, filters, bd_mode);
        end
    else
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
            filters = getDBfilter(vm, 0);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:2), 'db') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(3:end-1));
            filters = liftingfactortho(vm, 0, 1);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && ~strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end));
            filters = getDBfilter(vm, 1);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        elseif (strcmpi(wave_name(1:3), 'sym') && strcmpi(wave_name(end), 'x'))
            vm = str2double(wave_name(4:end-1));
            filters = liftingfactortho(vm, 1, 1);
            f = @(x, bd_mode) idwt_kernel_ortho(x, filters, bd_mode);
        end
    end
end

function x = dwt_kernel_97_dual(x,bd_mode)
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

function x = dwt_kernel_97(x,bd_mode)
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
    
function x = idwt_kernel_97_dual(x, bd_mode)
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
    
function x = idwt_kernel_97(x, bd_mode)
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
    
function x = dwt_kernel_53_dual(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x(1:2:end, :) = x(1:2:end, :)/alpha;
    x(2:2:end, :) = x(2:2:end, :)/beta;
    x = lifting_even_symm(lambda2, x, bd_mode);
    x = lifting_odd_symm(lambda1, x, bd_mode);
end
    
function x = dwt_kernel_53(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
    x = lifting_odd_symm(-lambda2, x, bd_mode);
    x = lifting_even_symm(-lambda1, x, bd_mode);
end
    
function x = idwt_kernel_53_dual(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x = lifting_odd_symm(-lambda1, x, bd_mode);
    x = lifting_even_symm(-lambda2, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)*alpha;
    x(2:2:end, :) = x(2:2:end, :)*beta;
end
   
function x = idwt_kernel_53(x, bd_mode)
    lambda1 = -1;
    lambda2 = 0.125;
    alpha = 2;
    beta = 0.5;
  
    x = lifting_even_symm(lambda1, x, bd_mode);
    x = lifting_odd_symm(lambda2, x, bd_mode);
    x(1:2:end, :) = x(1:2:end, :)/alpha;
    x(2:2:end, :) = x(2:2:end, :)/beta;
end
    
function x = dwt_kernel_pwl0_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
end
    
function x = dwt_kernel_pwl0(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
end
    
function x = idwt_kernel_pwl0_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_even_symm(-0.5, x, bd_mode);
end
    
function x = idwt_kernel_pwl0(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_odd_symm(0.5, x, bd_mode);
end

function x = dwt_kernel_pwl2_dual(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(0.5, x, bd_mode);
    x = lifting_odd_symm(-0.25, x, bd_mode);
end
    
function x = dwt_kernel_pwl2(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(-0.5, x, bd_mode);
    x = lifting_even_symm(0.25, x, bd_mode);
end
    
function x = idwt_kernel_pwl2_dual(x, bd_mode)
    x = x*sqrt(2);
    x = lifting_odd_symm(0.25, x, bd_mode);
    x = lifting_even_symm(-0.5, x, bd_mode);
end

function x = idwt_kernel_pwl2(x, bd_mode)
    x = x/sqrt(2);
    x = lifting_even_symm(-0.25, x, bd_mode);
    x = lifting_odd_symm(0.5, x, bd_mode);
end
    
function x = dwt_kernel_haar(x, bd_mode)
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
    
function x = idwt_kernel_haar(x, bd_mode)
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
    
function x = dwt_kernel_ortho_dual(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
    if bd_mode == 3
        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
    end
    if bd_mode >= 2
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
    if bd_mode >= 2
        x(1:size(filters.AL,2)) = x(1:size(filters.AL,2)) + y1;
        x((N-size(filters.AR,2)+1):N) = x((N-size(filters.AR,2)+1):N) + y2;
    end
end





function x = dwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
  
    y1 = 0; y2 = 0;
    if bd_mode == 3
        x(1:size(filters.A_L_pre,1)) = filters.A_L_pre*x(1:size(filters.A_L_pre,1));
        x((end-size(filters.A_R_pre,1)+1):end) = filters.A_R_pre*x((end-size(filters.A_R_pre,1)+1):end);
    end
    if bd_mode >= 2
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
    if bd_mode >= 2
        x(1:size(filters.AL,2)) = x(1:size(filters.AL,2)) + y1;
        x((N-size(filters.AR,2)+1):N) = x((N-size(filters.AR,2)+1):N) + y2;
    end
end




    
function x = idwt_kernel_ortho_dual(x, filters, bd_mode)
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if bd_mode >= 2
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
    
    if bd_mode >= 2
        x(1:size(filters.AL,1)) = x(1:size(filters.AL,1)) + y1;
        x((N-size(filters.AR,1)+1):N) = x((N-size(filters.AR,1)+1):N) + y2;
        if bd_mode == 3
            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
        end
    end
end


function x = idwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if bd_mode >= 2
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

    if bd_mode >= 2
        x(1:size(filters.AL,1)) = x(1:size(filters.AL,1)) + y1;
        x((N-size(filters.AR,1)+1):N) = x((N-size(filters.AR,1)+1):N) + y2;
        if bd_mode == 3
            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
        end
    end
end
    
function x=lifting_even_symm(lambda, x, bd_mode)
    N = size(x, 1);
    if ~bd_mode
        assert(mod(N,2) == 0)
    end
    if bd_mode
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
    if ~bd_mode
        assert(mod(N,2) == 0)
    end
    x(2:2:(N-1), :) = x(2:2:(N-1), :) + lambda*(x(1:2:(N-2), :) + x(3:2:N, :)); % This saves one multiplication
    if mod(N,2)==0 % last is even
        if bd_mode
            x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
        else
            x(N, :) = lambda*(x(1, :) + x(N-1, :)) + x(N, :);
        end
    end
end
    
function x = lifting_even(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    if bd_mode==0
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif bd_mode>=2
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
end
  
function x = lifting_odd(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    assert(mod(N,2) == 0)
    x(2:2:(N-1), :) = lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if bd_mode==0
        x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
    elseif bd_mode>=2
        x(N, :) = x(N, :) + lambda2*x(N-1, :);
    end
end