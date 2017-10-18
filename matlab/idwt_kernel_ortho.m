function x=idwt_kernel_ortho(x, filters, bd_mode)
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if ( strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'bd_pre') )
       y1 = filters.AL*x(1:size(filters.AL,2), :);
       y2 = filters.AR*x((N-size(filters.AR,2)+1):N, :);
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
        x(1:size(filters.AL,1), :) = x(1:size(filters.AL,1), :) + y1;
        x((N-size(filters.AR,1)+1):N, :) = x((N-size(filters.AR,1)+1):N, :) + y2;
%        if strcmpi(bd_mode, 'bd_pre')
%            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
%            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
%        end
    end
end

