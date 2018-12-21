function [f, prefilter]=find_kernel(wav_props, dual_wav_props, forward, dual, transpose, prefilter_mode)
    % Function which returns the default kernel of the library for use with the wavelet with properties encapsulated by the given parameters.
    % The kernel can be passed to the the internal functions (i)dwt1_impl_internal, (i)dwt2_impl_internal, (i)dwt3_impl_internal. 
    % User-defined kernels can also be passed to these internal functions: They simply have to take the x and bd_mode parameters, and return the 
    % transformed vectors.
    %
    % wav_props: Object which encapsulates the wavelet
    % dual_wav_props: Object which encapsulates the dual wavelet
    % forward: Whether the forward transform should be used
    % dual: (optional). Whether the dual wavelet should be applied. Default is 0
    % transpose: (optional). Default is 0
    % prefilter_mode: (optional). Default is 'none'
    
    if (~exist('dual','var')) dual = 0; end
    if (~exist('transpose','var')) transpose = 0; end
    if (~exist('prefilter_mode','var')) prefilter_mode = 'none'; end
    
    prefilter = @(x, forward) x;
    if transpose
        forward = ~forward;
        dual = ~dual;
    end
    if dual 
        temp = wav_props; wav_props = dual_wav_props; dual_wav_props = temp;
    end
    
    if strcmpi(wav_props.wave_name, 'Haar')
        if forward
            f = @dwt_kernel_haar;
        else
            f = @idwt_kernel_haar;
        end
    else
       if forward
           [f, prefilter]= find_kernel_dwt_general(wav_props, dual_wav_props, prefilter_mode);
       else
           [f, prefilter]= find_kernel_idwt_general(wav_props, dual_wav_props, prefilter_mode);
       end
    end
end

function [f, prefilter] = find_kernel_dwt_general(wav_props, dual_wav_props, prefilter_mode)
    prefilter = @(x, forward) x;
    if strcmpi(wav_props.wave_name, 'spline53') || strcmpi(wav_props.wave_name, 'cdf97') || strcmpi(wav_props.wave_name, 'pwl0') || strcmpi(wav_props.wave_name, 'pwl2')
        f = @(x, bd_mode) dwt_kernel_biortho(x, bd_mode, dual_wav_props);
        if strcmpi(prefilter_mode,'bd_pre')
            prefilter = @(x, forward) precond_impl(x, forward, wav_props);
        end
    elseif (strcmpi(wav_props.wave_name(1:2), 'db'))
        N = str2double(wav_props.wave_name(3:end));
        if N == 1
            f = @dwt_kernel_haar;
        else
            f = @(x, bd_mode) dwt_kernel_ortho(x, bd_mode, dual_wav_props);
            if strcmpi(prefilter_mode,'bd_pre')
                prefilter = @(x, forward) precond_impl(x, forward, wav_props);
            end
        end
    elseif strcmpi(wav_props.wave_name(1:3), 'sym')
        N = str2double(wav_props.wave_name(4:end));
        if N == 1
            f = @dwt_kernel_haar;
        else
            f = @(x, bd_mode) dwt_kernel_ortho(x, bd_mode, dual_wav_props);
            if strcmpi(prefilter_mode,'bd_pre')
                prefilter = @(x, forward) precond_impl(x, forward, wav_props);
            end
        end
    elseif strcmpi(wav_props.wave_name(1:6), 'spline')
        N = str2double(wav_props.wave_name(7));
        Ntilde = str2double(wav_props.wave_name(9));
        f = @(x, bd_mode) dwt_kernel_filters(x, bd_mode, dual_wav_props);
        if strcmpi(prefilter_mode,'bd_pre')
            prefilter = @(x, forward) precond_impl(x, forward, wav_props);
        end
    end
end

function [f, prefilter] = find_kernel_idwt_general(wav_props, dual_wav_props, prefilter_mode)
    prefilter = @(x, forward) x;
    if strcmpi(wav_props.wave_name, 'spline53') || strcmpi(wav_props.wave_name, 'cdf97') || strcmpi(wav_props.wave_name, 'pwl0') || strcmpi(wav_props.wave_name, 'pwl2')
        f = @(x, bd_mode) idwt_kernel_biortho(x, bd_mode, wav_props);
        if strcmpi(prefilter_mode,'bd_pre')
            prefilter = @(x, forward) precond_impl(x, forward, wav_props);
        end
    elseif strcmpi(wav_props.wave_name(1:2), 'db')
        N = str2double(wav_props.wave_name(3:end));
        if N == 1
            f = @idwt_kernel_haar;
        else
            f = @(x, bd_mode) idwt_kernel_ortho(x, bd_mode, wav_props);
            if strcmpi(prefilter_mode,'bd_pre')
                prefilter = @(x, forward) precond_impl(x, forward, wav_props);
            end
        end
    elseif strcmpi(wav_props.wave_name(1:3), 'sym')
        N = str2double(wav_props.wave_name(4:end));
        if N == 1
            f = @idwt_kernel_haar;
        else
            f = @(x, bd_mode) idwt_kernel_ortho(x, bd_mode, wav_props);
            if strcmpi(prefilter_mode,'bd_pre')
                prefilter = @(x, forward) precond_impl(x, forward, wav_props);
            end
        end
    elseif strcmpi(wav_props.wave_name(1:6), 'spline')
        N = str2double(wav_props.wave_name(7));
        Ntilde = str2double(wav_props.wave_name(9));
        f = @(x, bd_mode) idwt_kernel_filters(x, bd_mode, wav_props);
        if strcmpi(prefilter_mode,'bd_pre')
            prefilter = @(x, forward) precond_impl(x, forward, wav_props);
        end
    end
end

function x=dwt_kernel_filters(x, bd_mode, dual_wav_props)
    if strcmpi(bd_mode, 'bd')
        y1 = dual_wav_props.A_L'*x(1:size(dual_wav_props.A_L,1), :);
        y2 = dual_wav_props.A_R'*x((end-size(dual_wav_props.A_R,1)+1):end, :);
    end
    x0 = filter_impl(dual_wav_props.g0, x, bd_mode); 
    x1 = filter_impl(dual_wav_props.g1, x, bd_mode);
    x(1:2:end,:) = x0(1:2:end,:);
    x(2:2:end,:) = x1(2:2:end,:);
    if strcmpi(bd_mode, 'bd')
        x(1:size(dual_wav_props.A_L,2), :) ...
            = x(1:size(dual_wav_props.A_L,2), :) + y1;
        x((end-size(dual_wav_props.A_R,2)+1):end, :) ...
            = x((end-size(dual_wav_props.A_R,2)+1):end, :) + y2;
    end
end
% End dwt_kernel_filters
    
function x=idwt_kernel_filters(x, bd_mode, wav_props)
    if strcmpi(bd_mode, 'bd')
        y1 = wav_props.A_L*x(1:size(wav_props.A_L,2), :);
        y2 = wav_props.A_R*x((end-size(wav_props.A_R,2)+1):end, :);
    end
    x0 = x; x0(2:2:end,:) = 0;
    x1 = x; x1(1:2:end,:) = 0;
    x0 = filter_impl(wav_props.g0, x0, bd_mode);
    x1 = filter_impl(wav_props.g1, x1, bd_mode);
    x = x0 + x1;
    if strcmpi(bd_mode, 'bd')
        x(1:size(wav_props.A_L,1), :) ...
            = x(1:size(wav_props.A_L,1), :) + y1;
        x((end-size(wav_props.A_R,1)+1):end, :) ....
            = x((end-size(wav_props.A_R,1)+1):end, :) + y2;
    end
end
% End idwt_kernel_filters
    
function x=dwt_kernel_haar(x, bd_mode)
    x = x/sqrt(2);
    N = size(x, 1);
    if mod(N,2) == 1
        x(1:2, :) = [x(1, :) + x(2, :) - x(N, :);x(1, :) - x(2, :) - x(N, :)];
        x(N, :) = 2*x(N, :);
    else
        x(1:2, :) = [x(1, :) + x(2, :); x(1, :) - x(2, :)];
    end
    for k = 3:2:(N-1)
        x(k:(k+1), :) = [x(k, :) + x(k+1, :); x(k, :) - x(k+1, :)];
    end
end
% End dwt_kernel_haar

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
% End idwt_kernel_haar

function x=dwt_kernel_biortho(x, bd_mode, dual_wav_props)
    if strcmpi(bd_mode, 'bd')
        y1 = dual_wav_props.A_L'*x(1:size(dual_wav_props.A_L,1), :);
        y2 = dual_wav_props.A_R'*x((end-size(dual_wav_props.A_R,1)+1):end, :);
    end
    x(1:2:end, :) = x(1:2:end, :)/dual_wav_props.alpha;
    x(2:2:end, :) = x(2:2:end, :)/dual_wav_props.beta;
    iseven = ~dual_wav_props.last_even;
    for stepnr = (size(dual_wav_props.lambdas, 1)):(-1):1
        if iseven
            x = lifting_even_symm(dual_wav_props.lambdas(stepnr), x, bd_mode);
        else
            x = lifting_odd_symm(dual_wav_props.lambdas(stepnr), x, bd_mode);
        end
        iseven = ~iseven;
    end
    if strcmpi(bd_mode, 'bd')
        x(1:size(dual_wav_props.A_L,2), :) = x(1:size(dual_wav_props.A_L,2), :) + y1;
        x((end-size(dual_wav_props.A_R,2)+1):end, :) = x((end-size(dual_wav_props.A_R,2)+1):end, :) + y2;
    end
end

function x=idwt_kernel_biortho(x, bd_mode, wav_props)
    if strcmpi(bd_mode, 'bd')
        y1 = wav_props.A_L*x(1:size(wav_props.A_L,2), :);
        y2 = wav_props.A_R*x((end-size(wav_props.A_R,2)+1):end, :);
    end
    iseven = (mod(size(wav_props.lambdas, 1), 2) == wav_props.last_even);
    for stepnr = 1:(size(wav_props.lambdas, 1))
        if iseven
            x = lifting_even_symm(wav_props.lambdas(stepnr), x, bd_mode);
        else
            x = lifting_odd_symm(wav_props.lambdas(stepnr), x, bd_mode);
        end
        iseven = ~iseven;
    end
    x(1:2:end, :)=x(1:2:end, :)/wav_props.alpha;
    x(2:2:end, :)=x(2:2:end, :)/wav_props.beta;

    if strcmpi(bd_mode, 'bd')
        x(1:size(wav_props.A_L,1), :) = x(1:size(wav_props.A_L,1), :) + y1;
        x((end-size(wav_props.A_R,1)+1):end, :) = x((end-size(wav_props.A_R,1)+1):end, :) + y2;
    end
end

function x=dwt_kernel_ortho(x, bd_mode, dual_wav_props)
    if strcmpi(bd_mode, 'bd')
        y1 = dual_wav_props.A_L'*x(1:size(dual_wav_props.A_L,1), :);
        y2 = dual_wav_props.A_R'*x((end-size(dual_wav_props.A_R,1)+1):end, :);
    end
    x(1:2:end, :) = x(1:2:end, :)/dual_wav_props.alpha;
    x(2:2:end, :) = x(2:2:end, :)/dual_wav_props.beta;
    iseven = ~dual_wav_props.last_even;
    for stepnr = (size(dual_wav_props.lambdas, 1)):(-1):1
        if iseven
            x = lifting_even(dual_wav_props.lambdas(stepnr,2), ...
                             dual_wav_props.lambdas(stepnr,1), x, bd_mode);
        else
            x = lifting_odd(dual_wav_props.lambdas(stepnr,2), ...
                            dual_wav_props.lambdas(stepnr,1), x, bd_mode);
        end
        iseven = ~iseven;
    end
    if strcmpi(bd_mode, 'bd')
        x(1:size(dual_wav_props.A_L,2), :) = ...
            x(1:size(dual_wav_props.A_L,2), :) + y1;
        x((end-size(dual_wav_props.A_R,2)+1):end, :) = ...
            x((end-size(dual_wav_props.A_R,2)+1):end, :) + y2;
    end
end
% End dwt_kernel_ortho

function x=idwt_kernel_ortho(x, bd_mode, wav_props)    
    if strcmpi(bd_mode, 'bd')
       y1 = wav_props.A_L*x(1:size(wav_props.A_L,2), :);
       y2 = wav_props.A_R*x((end-size(wav_props.A_R,2)+1):end, :);
    end
    iseven = ( mod(size(wav_props.lambdas, 1), 2) == wav_props.last_even );
    for stepnr = 1:(size(wav_props.lambdas, 1))
        if iseven
            x = lifting_even(wav_props.lambdas(stepnr, 1), ...
                             wav_props.lambdas(stepnr, 2), x, bd_mode);
        else
            x = lifting_odd(wav_props.lambdas(stepnr, 1), ...
                            wav_props.lambdas(stepnr, 2), x, bd_mode);
        end
        iseven = ~iseven;
    end
    x(1:2:end, :)=x(1:2:end, :)/wav_props.alpha;
    x(2:2:end, :)=x(2:2:end, :)/wav_props.beta;

    if strcmpi(bd_mode, 'bd')
        x(1:size(wav_props.A_L,1), :) = ...
            x(1:size(wav_props.A_L,1), :) + y1;
        x((end-size(wav_props.A_R,1)+1):end, :) = ...
            x((end-size(wav_props.A_R,1)+1):end, :) + y2;
    end
end
% End idwt_kernel_ortho

function x=precond_impl(x, forward, wav_props)
    n = size(wav_props.A_L_pre_inv,1);
    if forward == 1
        x(1:n,:)           = wav_props.A_L_pre_inv\x(1:n,:);
        x((end-n+1):end,:) = wav_props.A_R_pre_inv\x((end-n+1):end,:);
    else
        x(1:n,:)           = wav_props.A_L_pre_inv*x(1:n,:);
        x((end-n+1):end,:) = wav_props.A_R_pre_inv*x((end-n+1):end,:);
    end
end

function x=lifting_even_symm(lambda, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        assert(mod(N,2) == 0)
    end
    if strcmpi(bd_mode, 'symm')
        x(1, :) = x(1, :) + 2*lambda*x(2, :); % Symmetric extension
    elseif strcmpi(bd_mode, 'per')
        x(1, :) = lambda*(x(2, :) + x(N, :)) + x(1, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        x(1, :) = lambda*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = x(3:2:(N-1),:) + lambda*(x(2:2:(N-2),:) + x(4:2:N,:));
    if mod(N,2) == 1 % last must also be included
        if strcmpi(bd_mode, 'symm')
            x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = x(N, :) + lambda*x(N-1, :);
        end
    end
end
% End lifting_even_symm
    
function x=lifting_odd_symm(lambda, x, bd_mode)
    N = size(x, 1);
    x(2:2:(N-1), :) = x(2:2:(N-1),:) + lambda*(x(1:2:(N-2),:) + x(3:2:N,:));
    if mod(N,2)==0 % last must also be included
        if strcmpi(bd_mode, 'symm')
            x(N, :) = x(N, :) + 2*lambda*x(N-1, :); % Symmetric extension
        elseif strcmpi(bd_mode, 'per')
            x(N, :) = lambda*(x(1, :) + x(N-1, :)) + x(N, :);
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = lambda*x(N-1, :) + x(N, :);
        end
    end
end
% End lifting_odd_symm

function x=lifting_even(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = ...
        lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        if mod(N,2) == 1 % last must also be included
            x(N, :) = x(N, :) + lambda2*x(N-1, :);
        end
    end
end
% End lifting_even

function x=lifting_odd(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    x(2:2:(N-1), :) = ...
        lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if mod(N,2) == 0 % last must also be included
        if strcmpi(bd_mode, 'per')
            x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = x(N, :) + lambda2*x(N-1, :);
        end
    end
end
% End lifting_odd